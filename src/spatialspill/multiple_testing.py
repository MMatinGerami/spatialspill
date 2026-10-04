"""Benjamini-Hochberg FDR control (thin wrapper, kept here to fix the method project-wide)."""

from __future__ import annotations

import numpy as np
from statsmodels.stats.multitest import multipletests


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaNs are passed through."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    m = np.isfinite(p)
    if m.any():
        out[m] = multipletests(p[m], method="fdr_bh")[1]
    return out
