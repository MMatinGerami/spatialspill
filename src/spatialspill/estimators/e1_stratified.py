"""E1: stratified difference in means with a stratified permutation null and BH-FDR.

For each target k, kind (autonomous or ring b) and recipient cell type c, the estimate is
the weighted average over strata s in c (strata = sample x cell type) of the mean difference
between treated and control cells, weights w_s = n_T n_C / (n_T + n_C). Standard errors use
within-group variances; p-values come from re-randomising guide labels within strata and
recomputing exposures and the statistic (so the null respects geometry and composition).
Two p-values are reported: the exact permutation p-value with the +1 correction
(``pvalue_perm``, floor 1/(n_perm+1)) and, by default, a permutation-calibrated z p-value
(``pvalue``): z = (observed - null mean) / null sd with a two-sided normal tail. The z version
has no floor and is what BH-FDR over tens of thousands of tests needs; its calibration is
checked on NTC pseudo-targets (NOTEBOOK amendment A1).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy import stats

from spatialspill.estimators.base import EstimateTable, Estimator
from spatialspill.estimators.groups import GroupConfig, group_masks
from spatialspill.exposure import Exposure
from spatialspill.permutation import permute_within_strata, strata_codes


def _stratum_indicator(strata: np.ndarray) -> sp.csr_matrix:
    S = strata.max() + 1
    return sp.csr_matrix(
        (np.ones(len(strata)), (np.arange(len(strata)), strata)), shape=(len(strata), S)
    )


def stratified_dim(
    T: sp.csr_matrix,
    C: sp.csr_matrix,
    Y: np.ndarray,
    strata: np.ndarray,
    stratum_group: np.ndarray,
    n_groups: int,
    min_cells: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Stratified difference in means for all targets and outcomes at once.

    Returns ``est, var, n_t, n_c`` with shapes (n_groups, K, G), (n_groups, K, G),
    (n_groups, K), (n_groups, K). ``stratum_group[s]`` maps each stratum to a reporting group
    (recipient cell type); group index ``n_groups - 1`` is "all" and pools every stratum.
    """
    S = _stratum_indicator(strata)  # n x S
    K, G = T.shape[1], Y.shape[1]
    nS = S.shape[1]
    Y2 = Y * Y
    # per-stratum, per-target sums: (S x K x G) via (K x n) @ (n_s x G) per stratum is costly;
    # use TS = T.multiply(S[:, s]) trick: loop over strata (S is small, tens to hundreds).
    est = np.zeros((n_groups, K, G))
    var = np.zeros((n_groups, K, G))
    wsum = np.zeros((n_groups, K))
    n_t = np.zeros((n_groups, K))
    n_c = np.zeros((n_groups, K))
    Tt = T.T.tocsr()
    Ct = C.T.tocsr()
    S_csc = S.tocsc()
    for s in range(nS):
        idx = S_csc.indices[S_csc.indptr[s] : S_csc.indptr[s + 1]]
        if len(idx) == 0:
            continue
        Ts = Tt[:, idx]
        Cs = Ct[:, idx]
        nt = np.asarray(Ts.sum(axis=1)).ravel()
        nc = np.asarray(Cs.sum(axis=1)).ravel()
        ok = (nt >= min_cells) & (nc >= min_cells)
        if not ok.any():
            continue
        Ys = Y[idx]
        Ys2 = Y2[idx]
        sum_t = np.asarray(Ts @ Ys)
        sum_c = np.asarray(Cs @ Ys)
        sq_t = np.asarray(Ts @ Ys2)
        sq_c = np.asarray(Cs @ Ys2)
        with np.errstate(invalid="ignore", divide="ignore"):
            mt = sum_t / nt[:, None]
            mc = sum_c / nc[:, None]
            vt = np.maximum(sq_t / nt[:, None] - mt**2, 0) * (nt / np.maximum(nt - 1, 1))[:, None]
            vc = np.maximum(sq_c / nc[:, None] - mc**2, 0) * (nc / np.maximum(nc - 1, 1))[:, None]
            d = mt - mc
            v = vt / nt[:, None] + vc / nc[:, None]
            w = nt * nc / (nt + nc)
        w = np.where(ok, w, 0.0)
        for grp in (stratum_group[s], n_groups - 1):
            est[grp] += w[:, None] * np.nan_to_num(d)
            var[grp] += (w[:, None] ** 2) * np.nan_to_num(v)
            wsum[grp] += w
            n_t[grp] += np.where(ok, nt, 0)
            n_c[grp] += np.where(ok, nc, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        est = est / wsum[:, :, None]
        var = var / (wsum[:, :, None] ** 2)
    return est, var, n_t, n_c


class E1Stratified(Estimator):
    name = "E1_stratified_dim"

    def __init__(
        self,
        n_perm: int = 200,
        seed: int = 0,
        strata_keys: tuple[str, ...] = ("sample", "cell_type"),
        min_cells: int = 5,
        groups: GroupConfig | None = None,
        report_by_cell_type: bool = True,
        pvalue: str = "z",
    ) -> None:
        if pvalue not in ("z", "perm"):
            raise ValueError("pvalue must be 'z' or 'perm'")
        self.pvalue = pvalue
        self.n_perm = n_perm
        self.seed = seed
        self.strata_keys = strata_keys
        self.min_cells = min_cells
        self.groups = groups or GroupConfig()
        self.report_by_cell_type = report_by_cell_type

    def _masks(self, exp: Exposure, is_ntc: np.ndarray, is_pert: np.ndarray) -> list[sp.csr_matrix]:
        T_auto, rings, C = group_masks(exp, is_ntc, is_pert, self.groups)
        return [T_auto, *rings, C]

    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        obs = adata.obs
        strata = strata_codes(obs, self.strata_keys)
        Y = np.asarray(outcomes, dtype=np.float64)
        is_ntc = obs["is_ntc"].to_numpy()
        is_pert = obs["is_perturbed"].to_numpy()
        # reporting groups: recipient cell type (+ "all")
        if self.report_by_cell_type and "cell_type" in self.strata_keys:
            ct_codes, ct_names = pd.factorize(obs["cell_type"].astype(str))
            # map stratum -> cell type via any member cell
            first = np.zeros(strata.max() + 1, dtype=int)
            first[strata] = np.arange(len(strata))
            stratum_group = ct_codes[first]
            groups = [*ct_names.tolist(), "all"]
        else:
            stratum_group = np.zeros(strata.max() + 1, dtype=int)
            groups = ["all"]
        n_groups = len(groups)

        masks = self._masks(exposure, is_ntc, is_pert)
        C = masks[-1]
        kinds = ["autonomous", *[f"ring{b}" for b in range(exposure.n_bins)]]
        obs_est, obs_var, obs_nt, obs_nc = [], [], [], []
        for M in masks[:-1]:
            e, v, nt, nc = stratified_dim(M, C, Y, strata, stratum_group, n_groups, self.min_cells)
            obs_est.append(e)
            obs_var.append(v)
            obs_nt.append(nt)
            obs_nc.append(nc)
        E = np.stack(obs_est)  # kinds x groups x K x G
        V = np.stack(obs_var)
        NT = np.stack(obs_nt)
        NC = np.stack(obs_nc)

        # permutation null: re-randomise labels within strata, recompute exposure + statistic
        exceed = np.zeros_like(E)
        n_valid = np.zeros_like(E)
        null_sum = np.zeros_like(E)
        null_sq = np.zeros_like(E)
        rng = np.random.default_rng(self.seed)
        labels0 = exposure.labels.copy()
        is_ntc_lab = np.array([t.startswith("NTC:") for t in labels0])
        for _ in range(self.n_perm):
            perm = permute_within_strata(labels0, strata, rng)
            # is_ntc / is_perturbed travel with the label
            p_ntc = np.array([t.startswith("NTC:") for t in perm])
            p_pert = np.array([exposure.is_perturbed_label.get(t, False) for t in perm])
            exposure.recompute(perm)
            pm = self._masks(exposure, p_ntc, p_pert)
            for ki, M in enumerate(pm[:-1]):
                e, _, _, _ = stratified_dim(
                    M, pm[-1], Y, strata, stratum_group, n_groups, self.min_cells
                )
                fin = np.isfinite(e) & np.isfinite(E[ki])
                exceed[ki] += fin & (np.abs(e) >= np.abs(E[ki]))
                n_valid[ki] += fin
                e0 = np.where(fin, e, 0.0)
                null_sum[ki] += e0
                null_sq[ki] += e0 * e0
        exposure.recompute(labels0)
        del is_ntc_lab
        pvals_perm = (exceed + 1) / (n_valid + 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            null_mean = null_sum / n_valid
            null_sd = np.sqrt(
                np.maximum(null_sq / n_valid - null_mean**2, 0)
                * n_valid
                / np.maximum(n_valid - 1, 1)
            )
            z = (E - null_mean) / null_sd
        pvals_z = 2 * stats.norm.sf(np.abs(z))
        pvals = pvals_z if self.pvalue == "z" else pvals_perm

        rows: list[dict[str, object]] = []
        se = np.sqrt(V)
        for ki, kind in enumerate(kinds):
            for gi, grp in enumerate(groups):
                for k, target in enumerate(exposure.targets):
                    ident = (
                        np.isfinite(E[ki, gi, k, 0])
                        and NT[ki, gi, k] >= self.min_cells
                        and NC[ki, gi, k] >= self.min_cells
                    )
                    for g, oname in enumerate(outcome_names):
                        est = E[ki, gi, k, g]
                        s = se[ki, gi, k, g]
                        rows.append(
                            {
                                "estimator": self.name,
                                "target": target,
                                "kind": "autonomous" if kind == "autonomous" else "spillover",
                                "ring": -1 if kind == "autonomous" else int(kind[4:]),
                                "cell_type": grp,
                                "outcome": oname,
                                "estimate": est,
                                "se": s,
                                "ci_low": est - 1.96 * s,
                                "ci_high": est + 1.96 * s,
                                "pvalue": pvals[ki, gi, k, g] if ident else np.nan,
                                "pvalue_perm": pvals_perm[ki, gi, k, g] if ident else np.nan,
                                "null_sd": null_sd[ki, gi, k, g] if ident else np.nan,
                                "n_treated": int(NT[ki, gi, k]),
                                "n_control": int(NC[ki, gi, k]),
                                "identified": bool(ident),
                            }
                        )
        tab = EstimateTable()
        tab.add(rows)
        return tab
