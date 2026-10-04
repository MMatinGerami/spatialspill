from __future__ import annotations

import re

import numpy as np
import pandas as pd

_NTC_RE = re.compile(r"(^|[_.-])(NTC|NT|non[-_ ]?target|safe|CTRL|control)", re.IGNORECASE)


def is_ntc_guide(name: str) -> bool:
    return bool(_NTC_RE.search(str(name)))


def guide_to_target(name: str) -> str:
    """Target gene symbol from heterogeneous guide names.

    Handles ``Gene.g2``, ``Gene_sg1``, ``sg54894_GENE``, ``sgrna12550_Gene``, ``NTC2``.
    """
    s = str(name).strip().strip('"')
    if is_ntc_guide(s):
        return "NTC"
    m = re.match(r"^sg(?:rna)?\d+[_-](.+)$", s, re.IGNORECASE)
    if m:
        return m.group(1)
    m = re.match(r"^(.+?)[._-](?:g|sg|sgRNA|guide)\d*$", s, re.IGNORECASE)
    if m:
        return m.group(1)
    return s


def dominant_guide(
    df: pd.DataFrame,
    unit_col: str,
    guide_col: str,
    count_col: str,
    min_count: int = 1,
    min_frac: float = 0.5,
) -> pd.DataFrame:
    """Per unit: dominant guide, its count, fraction of guide UMIs, number of distinct guides.

    Units whose dominant guide has fewer than ``min_count`` UMIs or less than ``min_frac`` of
    guide UMIs are labelled ``"none"`` (ambiguous or absent call). The fraction test is strict, so a
    two-guide tie at 0.5 is rejected.
    """
    g = df.groupby([unit_col, guide_col], observed=True)[count_col].sum().reset_index()
    tot = g.groupby(unit_col, observed=True)[count_col].sum().rename("total")
    nguides = g.groupby(unit_col, observed=True)[guide_col].nunique().rename("n_guides")
    top = g.sort_values(count_col, ascending=False).drop_duplicates(unit_col).set_index(unit_col)
    out = top.join(tot).join(nguides)
    out["frac"] = out[count_col] / out["total"]
    ok = (out[count_col] >= min_count) & (out["frac"] > min_frac)  # strict: ties are ambiguous
    out["guide"] = np.where(ok, out[guide_col].astype(str), "none")
    out["guide_confidence"] = out["frac"].where(ok, np.nan)
    return out[["guide", "guide_confidence", count_col, "total", "n_guides"]].rename(
        columns={count_col: "guide_umis", "total": "guide_umis_total"}
    )
