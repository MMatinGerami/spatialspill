"""Scoring estimator output against simulator truth: bias, interval coverage, power, AUROC."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import roc_auc_score


def truth_table(adata: AnnData, kind: str) -> pd.DataFrame:
    """Long table (target, outcome, true_lfc) for 'autonomous' or 'spillover' truth."""
    t = adata.uns["truth"]
    M = t["tau_auto"] if kind == "autonomous" else t["tau_spill"]
    targets = list(t["targets"])
    genes = list(adata.var_names)
    rows = [(tg, g, float(M[i, j])) for i, tg in enumerate(targets) for j, g in enumerate(genes)]
    return pd.DataFrame(rows, columns=["target", "outcome", "true_lfc"])


def score(
    est: pd.DataFrame, adata: AnnData, fdr: float = 0.1, cell_type: str = "all"
) -> pd.DataFrame:
    """Score one estimator table. Returns one row per (kind, ring).

    Bias and coverage are computed on the natural-log scale after converting the
    simulator's log fold change to the estimator's outcome scale by sign only: the
    estimators report differences of log1p-normalised means whose magnitude is not the
    planted LFC, so bias is reported on the sign-and-null structure (mean estimate for
    null pairs, which should be zero) and coverage as the fraction of null pairs whose
    95% interval contains zero. Power is the fraction of non-null pairs with q < fdr; AUROC
    ranks |estimate| / se against non-null status.
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
            null = m["true_lfc"] == 0
            z = (m["estimate"] / m["se"].replace(0, np.nan)).abs().fillna(0)
            auroc = roc_auc_score(~null, z) if 0 < (~null).sum() < len(m) else np.nan
            rows.append(
                {
                    "kind": kind,
                    "ring": int(s["ring"].iloc[0]),
                    "n_tests": len(m),
                    "n_nonnull": int((~null).sum()),
                    "null_mean_estimate": float(m.loc[null, "estimate"].mean()),
                    "null_coverage95": float(
                        ((m.loc[null, "ci_low"] <= 0) & (m.loc[null, "ci_high"] >= 0)).mean()
                    ),
                    "null_fpr_q": float((m.loc[null, "qvalue"] < fdr).mean()),
                    "null_fpr_p05": float((m.loc[null, "pvalue"] < 0.05).mean()),
                    "power_q": float((m.loc[~null, "qvalue"] < fdr).mean())
                    if (~null).any()
                    else np.nan,
                    "sign_agreement": float(
                        (
                            np.sign(m.loc[~null, "estimate"]) == np.sign(m.loc[~null, "true_lfc"])
                        ).mean()
                    )
                    if (~null).any()
                    else np.nan,
                    "auroc": float(auroc),
                    "fdp_q": float((m.loc[m["qvalue"] < fdr, "true_lfc"] == 0).mean())
                    if (m["qvalue"] < fdr).any()
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)
