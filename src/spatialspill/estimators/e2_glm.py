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

``n_perm`` > 0 (amendment A7) replaces the analytic t-test by the stratified permutation
null used in E1: labels are permuted within strata, exposures and regressions recomputed, and
the studentized coefficient is compared with its permutation distribution (permutation-
calibrated z p-value, exact permutation p in ``pvalue_perm``, CI rescaled by the null sd of t).
``cluster_tile_um`` > 0 clusters the robust standard errors by sample x spatial tile instead
of by sample, so that recipients sharing a niche are one cluster (amendment A6).
``spatial_basis`` > 0 adds, per sample, that many Gaussian radial basis functions centred on
k-means centres of the coordinates (bandwidth = median centre spacing unless given). This is
the pre-registered "spatial random effect" in spline form: it absorbs smooth niche variation
so that a clone's recipients are compared with controls at the same local level, without
splitting the data into tiles.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from spatialspill.estimators.base import EstimateTable, Estimator
from spatialspill.estimators.groups import GroupConfig, eligible_recipients
from spatialspill.exposure import Exposure
from spatialspill.permutation import permute_within_strata, strata_codes


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
        spatial_basis: int = 0,
        spatial_bandwidth_um: float | None = None,
        cluster_tile_um: float = 0.0,
        n_perm: int = 0,
        seed: int = 0,
        min_valid_perm: int | None = None,
    ) -> None:
        self.n_perm = n_perm
        self.seed = seed
        self.min_valid_perm = min_valid_perm if min_valid_perm is not None else max(20, n_perm // 2)
        self.cluster_tile_um = cluster_tile_um
        self.spatial_basis = spatial_basis
        self.spatial_bandwidth_um = spatial_bandwidth_um
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

    def _spatial_basis(self, adata: AnnData) -> np.ndarray:
        """Per-sample Gaussian radial basis on k-means centres (the spatial random effect)."""
        if self.spatial_basis <= 0:
            return np.zeros((adata.n_obs, 0))
        from sklearn.cluster import KMeans

        xy = np.asarray(adata.obsm["spatial"], dtype=float)
        samples = adata.obs["sample"].astype(str).to_numpy()
        blocks = []
        for s in np.unique(samples):
            idx = np.flatnonzero(samples == s)
            k = min(self.spatial_basis, max(1, len(idx) // 20))
            km = KMeans(n_clusters=k, n_init=1, random_state=0).fit(xy[idx])
            c = km.cluster_centers_
            d2 = ((xy[idx, None, :] - c[None, :, :]) ** 2).sum(-1)
            if self.spatial_bandwidth_um is None:
                cc = np.sort(((c[:, None, :] - c[None, :, :]) ** 2).sum(-1), axis=1)
                bw2 = float(np.median(cc[:, 1])) if k > 1 else float(d2.max())
            else:
                bw2 = self.spatial_bandwidth_um**2
            B = np.zeros((adata.n_obs, k))
            B[idx] = np.exp(-d2 / (2 * max(bw2, 1e-9)))
            blocks.append(B)
        return np.concatenate(blocks, axis=1)

    def _fit_one(
        self,
        idx: np.ndarray,
        own: np.ndarray,
        M: np.ndarray,
        strata: np.ndarray,
        W_all: np.ndarray,
        Y_all: np.ndarray,
        cluster_all: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, int, int, int, list[int]]:
        """OLS for one target and one selection; returns (B_coef, SE_coef, df, n_t, n_c, ring_n)
        where the coefficient rows are [own, ring_0, ..., ring_{B-1}]."""
        st = strata[idx]
        st_codes, st_inv = np.unique(st, return_inverse=True)
        S = np.zeros((len(idx), len(st_codes)))
        S[np.arange(len(idx)), st_inv] = 1.0
        W = W_all[idx]
        if W.shape[1]:
            means = (S.T @ W) / np.maximum(S.sum(0)[:, None], 1)
            W = W - S @ means
        X = np.column_stack([S, own[idx].astype(float), M[idx], W])
        n_t = int(own[idx].sum())
        n_c = int((~own[idx]).sum())
        ring_n = [int((M[idx][:, b] >= 1).sum()) for b in range(M.shape[1])]
        cluster: np.ndarray | None = cluster_all[idx] if self.cluster_by_sample else None
        if cluster is not None and len(np.unique(cluster)) < self.min_clusters:
            cluster = None
        B, SE, df = _ols_multi(X, Y_all[idx], cluster)
        p_own = S.shape[1]
        cols = [p_own, *range(p_own + 1, p_own + 1 + M.shape[1])]
        return B[cols], SE[cols], df, n_t, n_c, ring_n

    def _state(
        self, adata: AnnData, exposure: Exposure, labels: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[np.ndarray], np.ndarray]:
        """Label-dependent quantities for the current ``exposure`` state."""
        is_ntc = np.array([t.startswith("NTC:") for t in labels])
        is_pert = np.array([exposure.is_perturbed_label.get(t, False) for t in labels])
        elig = eligible_recipients(exposure, is_ntc, is_pert, self.groups)
        Z = exposure.label_matrix(labels).toarray() > 0
        tot = exposure.total_count_matrix().toarray()
        rings = [C.toarray() for C in exposure.counts]
        return elig, Z, tot, rings, exposure.any_perturbed_within.copy()

    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        obs = adata.obs
        Y_all = np.asarray(outcomes, dtype=np.float64)
        G = Y_all.shape[1]
        strata = strata_codes(obs, self.strata_keys)
        W_all = np.column_stack([self._covariates(adata, exposure), self._spatial_basis(adata)])
        is_pert_target = np.array(
            [exposure.is_perturbed_label.get(t, False) for t in exposure.targets]
        )
        ct = obs["cell_type"].astype(str).to_numpy()
        sample = obs["sample"].astype(str).to_numpy()
        if self.cluster_tile_um > 0:
            xy_all = np.asarray(adata.obsm["spatial"], dtype=float)
            tx = np.floor(xy_all[:, 0] / self.cluster_tile_um).astype(int)
            ty = np.floor(xy_all[:, 1] / self.cluster_tile_um).astype(int)
            cluster_all = np.array([f"{s_}|{a_}|{b_}" for s_, a_, b_ in zip(sample, tx, ty)])
        else:
            cluster_all = sample
        groups_out = [*sorted(set(ct)), "all"] if self.report_by_cell_type else ["all"]
        n_bins = exposure.n_bins
        n_coef = 1 + n_bins
        K, n_grp = len(exposure.targets), len(groups_out)

        def all_fits(
            labels: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            """Fit every (target, group); returns est, se, t, df, counts arrays
            with shapes (K, n_grp, n_coef, G) and counts (K, n_grp, n_coef, 2)."""
            elig, Z, tot, rings, any_pert = self._state(adata, exposure, labels)
            est = np.full((K, n_grp, n_coef, G), np.nan)
            se = np.full_like(est, np.nan)
            dfs = np.full((K, n_grp), np.nan)
            counts = np.zeros((K, n_grp, n_coef, 2), dtype=int)
            for k in range(K):
                own = Z[:, k]
                tk = tot[:, k]
                other = any_pert - (tk if is_pert_target[k] else 0.0)
                clean = (
                    other <= 0 if self.groups.clean_controls else np.ones(adata.n_obs, dtype=bool)
                )
                use = (own | (elig & ~own)) & clean
                if use.sum() < 2 * self.min_cells or own.sum() < self.min_cells:
                    continue
                M = np.column_stack([R[:, k] for R in rings])
                M[own] = 0.0  # ring coefficients from recipients only (A5)
                for gi, grp in enumerate(groups_out):
                    sel = use if grp == "all" else (use & (ct == grp))
                    if sel.sum() < 2 * self.min_cells:
                        continue
                    idx = np.flatnonzero(sel)
                    B, SE, df, n_t, n_c, ring_n = self._fit_one(
                        idx, own, M, strata, W_all, Y_all, cluster_all
                    )
                    est[k, gi] = B
                    se[k, gi] = SE
                    dfs[k, gi] = df
                    counts[k, gi, 0] = (n_t, n_c)
                    for b in range(n_bins):
                        counts[k, gi, 1 + b] = (ring_n[b], n_c - ring_n[b])
            with np.errstate(divide="ignore", invalid="ignore"):
                t = est / se
            return est, se, t, dfs, counts

        labels0 = exposure.labels.copy()
        est, se, t_obs, dfs, counts = all_fits(labels0)
        pvals = np.full_like(est, np.nan)
        null_sd = np.ones_like(est)
        enough = np.ones(est.shape[:3], dtype=bool)
        if self.n_perm > 0:
            rng = np.random.default_rng(self.seed)
            exceed = np.zeros_like(est)
            n_valid = np.zeros_like(est)
            s1 = np.zeros_like(est)
            s2 = np.zeros_like(est)
            for _ in range(self.n_perm):
                perm = permute_within_strata(labels0, strata, rng)
                exposure.recompute(perm)
                _, _, t_p, _, _ = all_fits(perm)
                fin = np.isfinite(t_p) & np.isfinite(t_obs)
                exceed += fin & (np.abs(t_p) >= np.abs(t_obs))
                n_valid += fin
                t0 = np.where(fin, t_p, 0.0)
                s1 += t0
                s2 += t0 * t0
            exposure.recompute(labels0)
            with np.errstate(invalid="ignore", divide="ignore"):
                mu = s1 / n_valid
                null_sd = np.sqrt(
                    np.maximum(s2 / n_valid - mu**2, 0) * n_valid / np.maximum(n_valid - 1, 1)
                )
                z = (t_obs - mu) / null_sd
            pvals = 2 * stats.norm.sf(np.abs(z))
            pvals_perm = (exceed + 1) / (n_valid + 1)
            enough = np.asarray((n_valid >= self.min_valid_perm).all(axis=3), dtype=bool)
            se = se * np.where(np.isfinite(null_sd), null_sd, 1.0)
            crit = np.full(dfs.shape, 1.96)
        else:
            for k in range(K):
                for gi in range(n_grp):
                    if np.isfinite(dfs[k, gi]):
                        pvals[k, gi] = 2 * stats.t.sf(np.abs(t_obs[k, gi]), dfs[k, gi])
            pvals_perm = np.full_like(est, np.nan)
            crit = np.where(
                np.isfinite(dfs), stats.t.ppf(0.975, np.where(np.isfinite(dfs), dfs, 1)), 1.96
            )

        rows: list[dict[str, object]] = []
        for k, target in enumerate(exposure.targets):
            for gi, grp in enumerate(groups_out):
                if not np.isfinite(dfs[k, gi]):
                    continue
                for ci in range(n_coef):
                    nt, nc = counts[k, gi, ci]
                    ident_base = (
                        nt >= self.min_cells and nc >= self.min_cells and bool(enough[k, gi, ci])
                    )
                    kind = "autonomous" if ci == 0 else "spillover"
                    ring = -1 if ci == 0 else ci - 1
                    for g, oname in enumerate(outcome_names):
                        e_ = float(est[k, gi, ci, g])
                        s_ = float(se[k, gi, ci, g])
                        pv = float(pvals[k, gi, ci, g])
                        ident = ident_base and np.isfinite(pv)
                        rows.append(
                            {
                                "estimator": self.name,
                                "target": target,
                                "kind": kind,
                                "ring": ring,
                                "cell_type": grp,
                                "outcome": oname,
                                "estimate": e_,
                                "se": s_,
                                "ci_low": e_ - float(crit[k, gi]) * s_,
                                "ci_high": e_ + float(crit[k, gi]) * s_,
                                "pvalue": pv if ident else np.nan,
                                "pvalue_perm": float(pvals_perm[k, gi, ci, g]) if ident else np.nan,
                                "n_treated": int(nt),
                                "n_control": int(nc),
                                "identified": bool(ident),
                            }
                        )
        tab = EstimateTable()
        tab.add(rows)
        return tab
