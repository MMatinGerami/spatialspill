"""Artifact detector (C3): how much apparent spillover is bleed-through, density or edge?

Bleed-through (docs/estimands.md, Section 6). If segmentation assigns a fraction alpha of cell
i's transcripts to a touching neighbour j, the first-ring "spillover" profile of target g
across genes is alpha times g's autonomous profile. For each target with an identified
autonomous profile a (vector over outcomes) and first-ring spillover profile s, the projection
coefficient alpha_hat = <s, a> / <a, a> and the explained fraction R2 = 1 - ||s - alpha a||^2 /
||s||^2 are computed; the residual r = s - alpha_hat a is the candidate biological spillover.
Outer rings serve as a control: bleed-through must vanish there, so alpha_hat(ring 0) -
alpha_hat(ring 2) is the per-target bleed-through estimate and the dataset-level artifact
fraction is the variance-weighted mean of R2 over targets in the first ring minus that in the
last ring. A positive alpha with R2 concentrated in ring 0 is the signature; biology that
mimics it (a knocked-out cell inducing the same programme in its neighbour) is not
distinguishable by expression alone and is listed as a limitation.

Density and edge effects are quantified by comparing estimates with and without the covariate
adjustment (E2 with and without density / edge distance); the function ``covariate_shift``
summarises the change.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _profiles(df: pd.DataFrame, kind: str, ring: int, cell_type: str) -> pd.DataFrame:
    sub = df[
        (df["kind"] == kind)
        & (df["ring"] == ring)
        & (df["cell_type"] == cell_type)
        & df["identified"]
    ]
    return sub.pivot_table(index="target", columns="outcome", values="estimate", aggfunc="first")


def bleedthrough_table(
    estimates: pd.DataFrame, cell_type: str = "all", min_outcomes: int = 10
) -> pd.DataFrame:
    """Per target and ring: alpha_hat, R2 and the number of outcomes used."""
    auto = _profiles(estimates, "autonomous", -1, cell_type)
    rings = sorted(estimates.loc[estimates["kind"] == "spillover", "ring"].unique())
    rows = []
    for b in rings:
        sp = _profiles(estimates, "spillover", int(b), cell_type)
        common_t = auto.index.intersection(sp.index)
        for t in common_t:
            cols = auto.columns.intersection(sp.columns)
            a = auto.loc[t, cols].to_numpy(dtype=float)
            s = sp.loc[t, cols].to_numpy(dtype=float)
            m = np.isfinite(a) & np.isfinite(s)
            if m.sum() < min_outcomes or np.sum(a[m] ** 2) == 0 or np.sum(s[m] ** 2) == 0:
                continue
            a, s = a[m], s[m]
            alpha = float(a @ s / (a @ a))
            resid = s - alpha * a
            r2 = float(1 - (resid @ resid) / (s @ s))
            rows.append(
                {
                    "target": t,
                    "ring": int(b),
                    "alpha_hat": alpha,
                    "r2": r2,
                    "n_outcomes": int(m.sum()),
                    "spill_norm": float(np.sqrt(s @ s)),
                    "auto_norm": float(np.sqrt(a @ a)),
                }
            )
    return pd.DataFrame(rows)


def artifact_fraction(bt: pd.DataFrame) -> dict[str, float | int | None]:
    """Dataset-level summary from :func:`bleedthrough_table`.

    ``r2_ring0_minus_last`` is the bleed-through artifact fraction: the share of first-ring
    spillover variance explained by the autonomous profile, beyond what the farthest ring shows.
    """
    if bt.empty:
        return {
            "n_targets": 0,
            "alpha_ring0": None,
            "alpha_last": None,
            "r2_ring0": None,
            "r2_last": None,
            "r2_ring0_minus_last": None,
        }
    rings = sorted(bt["ring"].unique())
    first, last = rings[0], rings[-1]

    def wmean(sub: pd.DataFrame, col: str) -> float | None:
        if sub.empty:
            return None
        w = sub["spill_norm"] ** 2
        return float(np.average(sub[col], weights=w)) if w.sum() > 0 else None

    b0 = bt[bt["ring"] == first]
    bl = bt[bt["ring"] == last]
    r0, rl = wmean(b0, "r2"), wmean(bl, "r2")
    return {
        "n_targets": int(bt["target"].nunique()),
        "alpha_ring0": wmean(b0, "alpha_hat"),
        "alpha_last": wmean(bl, "alpha_hat"),
        "r2_ring0": r0,
        "r2_last": rl,
        "r2_ring0_minus_last": (r0 - rl) if (r0 is not None and rl is not None) else None,
        "frac_targets_alpha_positive_ring0": float((b0["alpha_hat"] > 0).mean())
        if len(b0)
        else None,
    }


def covariate_shift(
    with_cov: pd.DataFrame, without_cov: pd.DataFrame, kind: str = "spillover"
) -> pd.DataFrame:
    """Compare estimates with and without covariate adjustment (same targets, rings, outcomes).

    Returns per ring: correlation of estimates, mean absolute change, and the fraction of
    hits (q < 0.1) in the unadjusted table that remain hits after adjustment.
    """
    keys = ["target", "ring", "cell_type", "outcome"]
    a = with_cov[with_cov["kind"] == kind][[*keys, "estimate", "qvalue", "identified"]]
    b = without_cov[without_cov["kind"] == kind][[*keys, "estimate", "qvalue", "identified"]]
    m = a.merge(b, on=keys, suffixes=("_adj", "_raw"))
    m = m[m["identified_adj"].astype(bool) & m["identified_raw"].astype(bool)]
    rows = []
    for ring, s in m.groupby("ring"):
        hits_raw = s["qvalue_raw"] < 0.1
        rows.append(
            {
                "ring": int(s["ring"].iloc[0]),
                "n": len(s),
                "corr": float(np.corrcoef(s["estimate_adj"], s["estimate_raw"])[0, 1])
                if len(s) > 2
                else np.nan,
                "mean_abs_change": float((s["estimate_adj"] - s["estimate_raw"]).abs().mean()),
                "hits_raw": int(hits_raw.sum()),
                "hits_retained": float((s.loc[hits_raw, "qvalue_adj"] < 0.1).mean())
                if hits_raw.any()
                else np.nan,
            }
        )
        del ring
    return pd.DataFrame(rows)
