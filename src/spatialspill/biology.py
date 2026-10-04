"""Ligand-receptor enrichment and niche-remodelling ranking (C6 support).

``load_lr_pairs`` reads OmniPath or CellPhoneDB tables into a set of (ligand, receptor) symbol
pairs. ``lr_enrichment`` asks whether inferred sender-gene to recipient-gene spillover pairs
(target g in the sender, outcome o in the recipient) are enriched for known axes where g is a
ligand (or receptor) and o is a cognate partner, by Fisher's exact test against all tested
pairs. ``rank_niche_genes`` aggregates a spillover table into one score per target.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def load_lr_pairs(
    omnipath_tsv: str | Path | None = None, cellphonedb_csv: str | Path | None = None
) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    if omnipath_tsv is not None:
        df = pd.read_csv(omnipath_tsv, sep="\t")
        cols = [
            c
            for c in df.columns
            if c.startswith("source_genesymbol") or c.startswith("target_genesymbol")
        ]
        if len(cols) >= 2:
            for a, b in zip(df[cols[0]].astype(str), df[cols[1]].astype(str)):
                for aa in a.split("_"):
                    for bb in b.split("_"):
                        pairs.add((aa.upper(), bb.upper()))
    if cellphonedb_csv is not None:
        df = pd.read_csv(cellphonedb_csv)
        ca = [c for c in df.columns if c in ("partner_a", "gene_name_a")]
        cb = [c for c in df.columns if c in ("partner_b", "gene_name_b")]
        if ca and cb:
            for a, b in zip(df[ca[0]].astype(str), df[cb[0]].astype(str)):
                pairs.add((a.upper(), b.upper()))
    return pairs


def lr_enrichment(
    spill: pd.DataFrame, lr: set[tuple[str, str]], q_threshold: float = 0.1
) -> dict[str, float | int]:
    """Fisher test: are significant (target, outcome) spillover pairs enriched for LR axes?"""
    s = spill[
        spill["identified"].astype(bool) & ~spill["target"].astype(str).str.startswith("NTC:")
    ].copy()
    s["t"] = s["target"].astype(str).str.upper()
    s["o"] = s["outcome"].astype(str).str.upper()
    s["is_lr"] = [((a, b) in lr) or ((b, a) in lr) for a, b in zip(s["t"], s["o"])]
    s["hit"] = s["qvalue"] < q_threshold
    a = int((s["hit"] & s["is_lr"]).sum())
    b = int((s["hit"] & ~s["is_lr"]).sum())
    c = int((~s["hit"] & s["is_lr"]).sum())
    d = int((~s["hit"] & ~s["is_lr"]).sum())
    if a + b == 0 or a + c == 0:
        return {
            "n_pairs": len(s),
            "n_lr_pairs": a + c,
            "n_hits": a + b,
            "hits_lr": a,
            "odds_ratio": np.nan,
            "pvalue": np.nan,
        }
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "n_pairs": len(s),
        "n_lr_pairs": a + c,
        "n_hits": a + b,
        "hits_lr": a,
        "odds_ratio": float(odds),
        "pvalue": float(p),
    }


def rank_niche_genes(
    spill: pd.DataFrame, autonomous: pd.DataFrame | None = None, cell_type: str = "all"
) -> pd.DataFrame:
    """One row per target: aggregate spillover evidence, with the autonomous profile for context.

    score = mean over identified (ring, outcome) of z^2 minus 1 (expected zero under the null),
    so a target with no spillover has score about 0; n_tests and the best q-value are kept.
    """
    s = spill[
        (spill["kind"] == "spillover")
        & (spill["cell_type"] == cell_type)
        & spill["identified"].astype(bool)
    ]
    s = s[~s["target"].astype(str).str.startswith("NTC:")].copy()
    s["z2"] = (s["estimate"] / s["se"].replace(0, np.nan)) ** 2
    g = s.groupby("target").agg(
        n_tests=("z2", "size"),
        score=("z2", lambda x: float(np.nanmean(x) - 1)),
        best_q=("qvalue", "min"),
        n_hits=("qvalue", lambda q: int((q < 0.1).sum())),
    )
    ntc = spill[
        (spill["kind"] == "spillover")
        & (spill["cell_type"] == cell_type)
        & spill["identified"].astype(bool)
        & spill["target"].astype(str).str.startswith("NTC:")
    ].copy()
    if len(ntc):
        ntc["z2"] = (ntc["estimate"] / ntc["se"].replace(0, np.nan)) ** 2
        ntc_scores = ntc.groupby("target")["z2"].apply(lambda x: float(np.nanmean(x) - 1))
        g["ntc_score_q95"] = float(np.nanquantile(ntc_scores, 0.95)) if len(ntc_scores) else np.nan
    if autonomous is not None:
        a = autonomous[
            (autonomous["kind"] == "autonomous")
            & (autonomous["cell_type"] == cell_type)
            & autonomous["identified"].astype(bool)
        ]
        g["autonomous_hits"] = a.groupby("target")["qvalue"].apply(lambda q: int((q < 0.1).sum()))
    return g.sort_values("score", ascending=False).reset_index()
