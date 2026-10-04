"""E2: regression adjustment with covariates and distance-bin exposures.

For each target g and each outcome, a linear model on the log1p-normalised outcome (fast,
robust; the negative-binomial GLM on counts is available as ``family="nb"`` for a subset of
outcomes) is fitted on the eligible cells (recipients/controls per ADR-003 plus the g cells):

    y_i = sum_s a_s 1{stratum_i = s} + b * own_i + sum_b c_b m_i^b + d' W_i + e_i

where ``own_i`` indicates the cell carries g (autonomous effect b), ``m_i^b`` is the number of
g cells in ring b (per-neighbour spillover c_b), and W holds local density, log area and edge
distance (centred within stratum). Standard errors are heteroskedasticity-robust (HC1) and,
when more than one sample is present, clustered by sample. P-values are from the t
distribution (cluster degrees of freedom when at least ``min_clusters`` samples are present,
else HC1 with residual degrees of freedom); the permutation null of E1 is not repeated here, so the mandatory NTC
calibration check for E2 is the empirical false-positive rate over NTC pseudo-targets.

Rings beyond D_max (far field) can be added through ``far_bins_um`` as a built-in negative
control: their coefficients must be null.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from spatialspill.estimators.base import EstimateTable, Estimator
from spatialspill.estimators.groups import GroupConfig, eligible_recipients
from spatialspill.exposure import Exposure
from spatialspill.permutation import strata_codes


def _ols_multi(
    X: np.ndarray, Y: np.ndarray, cluster: np.ndarray | None
) -> tuple[np.ndarray, np.ndarray, int]:
    """OLS of many outcomes on one design. Returns (beta, se, df) with HC1 or cluster-robust SE."""
    n, p = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ (X.T @ Y)  # p x G
    R = Y - X @ B  # n x G
    if cluster is None or len(np.unique(cluster)) < 2:
        # HC1: (X'X)^-1 X' diag(r^2) X (X'X)^-1 * n/(n-p), per outcome
        se = np.empty_like(B)
        for g in range(Y.shape[1]):
            r2 = R[:, g] ** 2
            meat = (X * r2[:, None]).T @ X
            V = XtX_inv @ meat @ XtX_inv * (n / max(n - p, 1))
            se[:, g] = np.sqrt(np.maximum(np.diag(V), 0))
        df = max(n - p, 1)
    else:
        groups, inv = np.unique(cluster, return_inverse=True)
        Gn = len(groups)
        se = np.empty_like(B)
        for g in range(Y.shape[1]):
            U = np.zeros((Gn, p))
            np.add.at(U, inv, X * R[:, g][:, None])
            meat = U.T @ U
            V = XtX_inv @ meat @ XtX_inv * (Gn / max(Gn - 1, 1)) * ((n - 1) / max(n - p, 1))
            se[:, g] = np.sqrt(np.maximum(np.diag(V), 0))
        df = max(Gn - 1, 1)
    return B, se, df


class E2GLM(Estimator):
    name = "E2_glm"

    def __init__(
        self,
        strata_keys: tuple[str, ...] = ("sample", "cell_type"),
        covariates: tuple[str, ...] = ("local_density", "log_area", "edge_distance"),
        min_cells: int = 5,
        groups: GroupConfig | None = None,
        cluster_by_sample: bool = True,
        min_clusters: int = 5,
        report_by_cell_type: bool = True,
    ) -> None:
        self.strata_keys = strata_keys
        self.covariates = covariates
        self.min_cells = min_cells
        self.groups = groups or GroupConfig()
        self.cluster_by_sample = cluster_by_sample
        self.min_clusters = min_clusters
        self.report_by_cell_type = report_by_cell_type

    def _covariates(self, adata: AnnData, exposure: Exposure) -> np.ndarray:
        obs = adata.obs
        cols = []
        for c in self.covariates:
            if c == "local_density":
                v = exposure.n_neighbors_within.astype(float)
            elif c == "log_area":
                a = pd.to_numeric(obs["area"], errors="coerce").to_numpy(dtype=float)
                v = np.log(np.where(np.isfinite(a) & (a > 0), a, np.nan))
                if np.isnan(v).all():
                    continue
                v = np.where(np.isnan(v), np.nanmean(v), v)
            elif c in obs:
                v = pd.to_numeric(obs[c], errors="coerce").fillna(0).to_numpy(dtype=float)
            else:
                continue
            cols.append(v)
        return np.column_stack(cols) if cols else np.zeros((adata.n_obs, 0))

    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        obs = adata.obs
        Y_all = np.asarray(outcomes, dtype=np.float64)
        strata = strata_codes(obs, self.strata_keys)
        is_ntc = obs["is_ntc"].to_numpy()
        is_pert = obs["is_perturbed"].to_numpy()
        elig = eligible_recipients(exposure, is_ntc, is_pert, self.groups)
        W_all = self._covariates(adata, exposure)
        Z = exposure.label_matrix(exposure.labels).tocsc()
        tot = exposure.total_count_matrix().tocsc()
        ring_counts = [C.tocsc() for C in exposure.counts]
        any_pert = exposure.any_perturbed_within
        is_pert_target = np.array(
            [exposure.is_perturbed_label.get(t, False) for t in exposure.targets]
        )
        ct = obs["cell_type"].astype(str).to_numpy()
        sample = obs["sample"].astype(str).to_numpy()
        groups_out = [*sorted(set(ct)), "all"] if self.report_by_cell_type else ["all"]
        rows: list[dict[str, object]] = []
        n_bins = exposure.n_bins

        for k, target in enumerate(exposure.targets):
            own = np.asarray(Z[:, k].todense()).ravel() > 0
            tk = np.asarray(tot[:, k].todense()).ravel()
            other = any_pert - (tk if is_pert_target[k] else 0.0)
            clean = other <= 0 if self.groups.clean_controls else np.ones(adata.n_obs, dtype=bool)
            # own-target cells are never recipients or controls for their own target
            use = (own | (elig & ~own)) & clean
            if use.sum() < 2 * self.min_cells or own.sum() < self.min_cells:
                continue
            M = np.column_stack(
                [np.asarray(C[:, k].todense()).ravel() for C in ring_counts]
            )  # n x B
            for grp in groups_out:
                sel = use if grp == "all" else (use & (ct == grp))
                if sel.sum() < 2 * self.min_cells:
                    continue
                idx = np.flatnonzero(sel)
                st = strata[idx]
                st_codes, st_inv = np.unique(st, return_inverse=True)
                S = np.zeros((len(idx), len(st_codes)))
                S[np.arange(len(idx)), st_inv] = 1.0
                W = W_all[idx]
                # centre covariates within stratum
                if W.shape[1]:
                    means = (S.T @ W) / np.maximum(S.sum(0)[:, None], 1)
                    W = W - S @ means
                X = np.column_stack([S, own[idx].astype(float), M[idx], W])
                n_t = int(own[idx].sum())
                n_c = int((~own[idx]).sum())
                ring_n = [int((M[idx][:, b] >= 1).sum()) for b in range(n_bins)]
                cluster = sample[idx] if self.cluster_by_sample else None
                if cluster is not None and len(np.unique(cluster)) < self.min_clusters:
                    cluster = None  # too few clusters for a cluster-robust SE; fall back to HC1
                B, SE, df = _ols_multi(X, Y_all[idx], cluster)
                p_own = S.shape[1]
                coefs = {"autonomous": (p_own, -1, n_t, n_c)}
                for b in range(n_bins):
                    coefs[f"ring{b}"] = (p_own + 1 + b, b, ring_n[b], n_c - ring_n[b])
                for kind, (col, ring, nt, nc) in coefs.items():
                    ident = nt >= self.min_cells and nc >= self.min_cells
                    est = B[col]
                    se = SE[col]
                    with np.errstate(divide="ignore", invalid="ignore"):
                        tstat = est / se
                    pv = 2 * stats.t.sf(np.abs(tstat), df)
                    for g, oname in enumerate(outcome_names):
                        rows.append(
                            {
                                "estimator": self.name,
                                "target": target,
                                "kind": "autonomous" if kind == "autonomous" else "spillover",
                                "ring": ring,
                                "cell_type": grp,
                                "outcome": oname,
                                "estimate": float(est[g]),
                                "se": float(se[g]),
                                "ci_low": float(est[g] - stats.t.ppf(0.975, df) * se[g]),
                                "ci_high": float(est[g] + stats.t.ppf(0.975, df) * se[g]),
                                "pvalue": float(pv[g]) if ident and np.isfinite(pv[g]) else np.nan,
                                "n_treated": nt,
                                "n_control": nc,
                                "identified": bool(ident and np.isfinite(pv[g])),
                            }
                        )
        tab = EstimateTable()
        tab.add(rows)
        return tab
