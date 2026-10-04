"""Scoring estimator output against simulator truth: bias, interval coverage, power, AUROC."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import roc_auc_score


def truth_table(adata: AnnData, kind: str, eps: float = 0.05) -> pd.DataFrame:
    """Long table (target, outcome, true_lfc, effective_lfc, status) for one effect kind.

    ``effective_lfc`` is the planted log fold change corrected for the compositional shift that
    size-factor normalisation induces: for a cell with baseline gene proportions p and planted
    LFC vector tau, the normalised log change of gene g is tau_g - log(sum_j p_j exp(tau_j)).
    ``status`` is "nonnull" (planted tau != 0), "null" (tau == 0 and |effective| < eps) or
    "compositional" (tau == 0 but |effective| >= eps; excluded from null scoring, counted).
    NTC pseudo-targets are exact nulls.
    """
    t = adata.uns["truth"]
    M = np.asarray(t["tau_auto"] if kind == "autonomous" else t["tau_spill"], dtype=float)
    p = np.asarray(t.get("base_props", np.ones(M.shape[1]) / M.shape[1]), dtype=float)
    targets = list(t["targets"])
    genes = list(adata.var_names)
    rows = []
    for i, tg in enumerate(targets):
        shift = np.log(np.sum(p * np.exp(M[i])))
        eff = M[i] - shift
        for j, g in enumerate(genes):
            if M[i, j] != 0:
                status = "nonnull"
            elif abs(eff[j]) < eps:
                status = "null"
            else:
                status = "compositional"
            rows.append((tg, g, float(M[i, j]), float(eff[j]), status))
    ntc_guides = sorted(set(adata.obs.loc[adata.obs["is_ntc"], "guide"].astype(str)))
    rows += [(f"NTC:{gd}", g, 0.0, 0.0, "null") for gd in ntc_guides for g in genes]
    return pd.DataFrame(rows, columns=["target", "outcome", "true_lfc", "effective_lfc", "status"])


def score(
    est: pd.DataFrame, adata: AnnData, fdr: float = 0.1, cell_type: str = "all"
) -> pd.DataFrame:
    """Score one estimator table. Returns one row per (kind, ring).

    Pairs are "null", "nonnull" or "compositional" (see :func:`truth_table`). Null
    calibration (mean estimate, 95% coverage of zero, FPR at q and at p < 0.05) uses null
    pairs; NTC pseudo-targets are scored separately as exact nulls. Power is the fraction of
    non-null pairs with q < fdr; sign agreement is against the composition-corrected
    effective LFC; AUROC ranks |estimate| / se against non-null status with compositional
    pairs excluded; FDP is the share of null pairs among hits, compositional pairs excluded.
    """
    rows = []
    e = est[est["cell_type"] == cell_type]
    for kind in ("autonomous", "spillover"):
        tt = truth_table(adata, kind)
        sub = e[e["kind"] == kind]
        if sub.empty:
            continue
        for _ring, s in sub.groupby("ring"):
            m = s.merge(tt, on=["target", "outcome"], how="inner")
            m = m[m["identified"].astype(bool)]
            if m.empty:
                continue
            null = m["status"] == "null"
            nonnull = m["status"] == "nonnull"
            comp = m["status"] == "compositional"
            ntc_m = m["target"].astype(str).str.startswith("NTC:")
            z = (m["estimate"] / m["se"].replace(0, np.nan)).abs().fillna(0)
            mm = m[~comp]
            zz = z[~comp]
            auroc = (
                roc_auc_score(nonnull[~comp], zz) if 0 < nonnull[~comp].sum() < len(mm) else np.nan
            )
            rows.append(
                {
                    "kind": kind,
                    "ring": int(s["ring"].iloc[0]),
                    "n_tests": len(m),
                    "n_nonnull": int(nonnull.sum()),
                    "n_compositional": int(comp.sum()),
                    "null_mean_estimate": float(m.loc[null, "estimate"].mean()),
                    "null_coverage95": float(
                        ((m.loc[null, "ci_low"] <= 0) & (m.loc[null, "ci_high"] >= 0)).mean()
                    ),
                    "null_fpr_q": float((m.loc[null, "qvalue"] < fdr).mean()),
                    "null_fpr_p05": float((m.loc[null, "pvalue"] < 0.05).mean()),
                    "ntc_coverage95": float(
                        ((m.loc[ntc_m, "ci_low"] <= 0) & (m.loc[ntc_m, "ci_high"] >= 0)).mean()
                    )
                    if ntc_m.any()
                    else np.nan,
                    "ntc_fpr_q": float((m.loc[ntc_m, "qvalue"] < fdr).mean())
                    if ntc_m.any()
                    else np.nan,
                    "ntc_fpr_p05": float((m.loc[ntc_m, "pvalue"] < 0.05).mean())
                    if ntc_m.any()
                    else np.nan,
                    "power_q": float((m.loc[nonnull, "qvalue"] < fdr).mean())
                    if nonnull.any()
                    else np.nan,
                    "sign_agreement": float(
                        (
                            np.sign(m.loc[nonnull, "estimate"])
                            == np.sign(m.loc[nonnull, "effective_lfc"])
                        ).mean()
                    )
                    if nonnull.any()
                    else np.nan,
                    "auroc": float(auroc),
                    "fdp_q": float((m.loc[(m["qvalue"] < fdr) & ~comp, "status"] == "null").mean())
                    if ((m["qvalue"] < fdr) & ~comp).any()
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)
