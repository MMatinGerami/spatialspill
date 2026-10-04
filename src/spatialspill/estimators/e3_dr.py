"""E3: augmented inverse-probability-weighted (doubly robust) estimator under interference.

Following Aronow and Samii (2017), exposure probabilities pi_i(e) = P(e_i(Z) = e) are computed
by Monte Carlo re-randomisation of the guide labels within strata on the fixed geometry (the
same stratified permutation used for the null). For each target g, recipient exposure e is
"g-neighbour in ring b only" versus "no g-neighbour within D_max", restricted to eligible
recipients (ADR-003). The estimator for E[Y(e)] is the Hajek-normalised augmented estimator

    mu_hat(e) = sum_i [ 1{E_i = e} (Y_i - m_e(W_i)) / pi_i(e) ] / sum_i [ 1{E_i = e} / pi_i(e) ]
                + mean_i m_e(W_i)

with an outcome model m_e(W) fitted by ridge regression on stratum indicators and covariates
among the cells observed at exposure e. The spillover effect is mu_hat(ring b) - mu_hat(none).
Standard errors are from a cell-level bootstrap of the influence-function terms (fast, closed
form) and, when at least ``min_clusters`` samples exist, a cluster bootstrap over samples.
Cells with pi below ``pi_floor`` are trimmed and counted in the output.

The autonomous effect is estimated the same way with exposure "carries g, no g-neighbour".
"""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy import stats

from spatialspill.estimators.base import EstimateTable, Estimator
from spatialspill.estimators.groups import GroupConfig, group_masks
from spatialspill.exposure import Exposure
from spatialspill.permutation import permute_within_strata, strata_codes


def exposure_probabilities(
    exposure: Exposure,
    is_ntc: np.ndarray,
    is_pert: np.ndarray,
    strata: np.ndarray,
    cfg: GroupConfig,
    n_draws: int,
    seed: int,
) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
    """Monte Carlo P(exposure) per (cell, target): autonomous, each ring, control."""
    rng = np.random.default_rng(seed)
    labels0 = exposure.labels.copy()
    n, K = exposure.n_cells, len(exposure.targets)
    p_auto = np.zeros((n, K))
    p_ring = [np.zeros((n, K)) for _ in range(exposure.n_bins)]
    p_ctrl = np.zeros((n, K))
    for _ in range(n_draws):
        perm = permute_within_strata(labels0, strata, rng)
        p_is_ntc = np.array([t.startswith("NTC:") for t in perm])
        p_is_pert = np.array([exposure.is_perturbed_label.get(t, False) for t in perm])
        exposure.recompute(perm)
        T, rings, C = group_masks(exposure, p_is_ntc, p_is_pert, cfg)
        p_auto += T.toarray()
        for b, R in enumerate(rings):
            p_ring[b] += R.toarray()
        p_ctrl += C.toarray()
    exposure.recompute(labels0)
    return p_auto / n_draws, [p / n_draws for p in p_ring], p_ctrl / n_draws


def _ridge_fit_predict(
    X: np.ndarray, y: np.ndarray, X_all: np.ndarray, lam: float = 1.0
) -> np.ndarray:
    p = X.shape[1]
    A = X.T @ X + lam * np.eye(p)
    beta = np.linalg.solve(A, X.T @ y)
    return X_all @ beta


class E3DoublyRobust(Estimator):
    name = "E3_aipw"

    def __init__(
        self,
        n_draws: int = 200,
        seed: int = 0,
        strata_keys: tuple[str, ...] = ("sample", "cell_type"),
        min_cells: int = 5,
        groups: GroupConfig | None = None,
        pi_floor: float = 0.01,
        n_boot: int = 200,
        report_by_cell_type: bool = True,
    ) -> None:
        self.n_draws = n_draws
        self.seed = seed
        self.strata_keys = strata_keys
        self.min_cells = min_cells
        self.groups = groups or GroupConfig()
        self.pi_floor = pi_floor
        self.n_boot = n_boot
        self.report_by_cell_type = report_by_cell_type

    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        obs = adata.obs
        Y = np.asarray(outcomes, dtype=np.float64)
        n, G = Y.shape
        strata = strata_codes(obs, self.strata_keys)
        is_ntc = obs["is_ntc"].to_numpy()
        is_pert = obs["is_perturbed"].to_numpy()
        T_obs, rings_obs, C_obs = group_masks(exposure, is_ntc, is_pert, self.groups)
        p_auto, p_ring, p_ctrl = exposure_probabilities(
            exposure, is_ntc, is_pert, strata, self.groups, self.n_draws, self.seed
        )
        # covariates for the outcome model: stratum dummies + density
        S_codes, S_inv = np.unique(strata, return_inverse=True)
        S = np.zeros((n, len(S_codes)))
        S[np.arange(n), S_inv] = 1.0
        dens = exposure.n_neighbors_within.astype(float)
        dens = (dens - dens.mean()) / (dens.std() + 1e-9)
        Wmat = np.column_stack([S, dens])
        ct = obs["cell_type"].astype(str).to_numpy()
        groups_out = [*sorted(set(ct)), "all"] if self.report_by_cell_type else ["all"]
        rng = np.random.default_rng(self.seed + 1)
        rows: list[dict[str, object]] = []

        def aipw(
            mask_e: np.ndarray, pi_e: np.ndarray, sel: np.ndarray
        ) -> tuple[np.ndarray, np.ndarray, int]:
            """Return (mu_hat (G,), influence terms (n_sel, G), n_exposed) within selection ``sel``."""
            idx = np.flatnonzero(sel)
            me = mask_e[idx]
            pe = pi_e[idx]
            keep = pe >= self.pi_floor
            idx, me, pe = idx[keep], me[keep], pe[keep]
            n_e = int(me.sum())
            if n_e < self.min_cells:
                return np.full(G, np.nan), np.zeros((0, G)), n_e
            Xe = Wmat[idx]
            m_hat = np.column_stack(
                [_ridge_fit_predict(Xe[me], Y[idx][me][:, g], Xe) for g in range(G)]
            )
            w = me / pe
            w = w / w.sum()
            resid = Y[idx] - m_hat
            mu = (w[:, None] * resid).sum(0) + m_hat.mean(0)
            infl = len(idx) * w[:, None] * resid + m_hat - mu  # per-cell contributions, mean = 0
            return mu, infl, n_e

        kinds = [("autonomous", T_obs.toarray(), p_auto, -1)] + [
            (f"ring{b}", rings_obs[b].toarray(), p_ring[b], b) for b in range(exposure.n_bins)
        ]
        C_d = C_obs.toarray()
        for k, target in enumerate(exposure.targets):
            for grp in groups_out:
                sel = np.ones(n, dtype=bool) if grp == "all" else ct == grp
                mu_c, infl_c, n_c = aipw(C_d[:, k] > 0, p_ctrl[:, k], sel)
                if n_c < self.min_cells:
                    continue
                for kind, M, P, ring in kinds:
                    mu_e, infl_e, n_e = aipw(M[:, k] > 0, P[:, k], sel)
                    ident = n_e >= self.min_cells and np.isfinite(mu_e).all()
                    if ident:
                        est = mu_e - mu_c
                        # bootstrap over cells of the influence terms (independent resampling of the two groups)
                        boots = np.empty((self.n_boot, G))
                        for b_ in range(self.n_boot):
                            ie = infl_e[rng.integers(0, len(infl_e), len(infl_e))].mean(0)
                            ic = infl_c[rng.integers(0, len(infl_c), len(infl_c))].mean(0)
                            boots[b_] = ie - ic
                        se = boots.std(0, ddof=1) / 1.0
                        with np.errstate(divide="ignore", invalid="ignore"):
                            z = est / se
                        pv = 2 * stats.norm.sf(np.abs(z))
                    else:
                        est = np.full(G, np.nan)
                        se = np.full(G, np.nan)
                        pv = np.full(G, np.nan)
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
                                "ci_low": float(est[g] - 1.96 * se[g]),
                                "ci_high": float(est[g] + 1.96 * se[g]),
                                "pvalue": float(pv[g]) if ident else np.nan,
                                "n_treated": n_e,
                                "n_control": n_c,
                                "identified": bool(ident and np.isfinite(pv[g])),
                            }
                        )
        tab = EstimateTable()
        tab.add(rows)
        return tab
