"""Outcome matrices on the analysis scale."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from anndata import AnnData


def lognorm(
    adata: AnnData, target_sum: float | None = None, layer: str | None = None
) -> np.ndarray:
    """log1p of size-factor-normalised counts, dense float32.

    ``target_sum`` defaults to the median total count, as in scanpy.
    """
    X = adata.layers[layer] if layer else adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = X.astype(np.float64)
    tot = X.sum(axis=1)
    if target_sum is None:
        target_sum = float(np.median(tot[tot > 0])) if (tot > 0).any() else 1.0
    sf = np.where(tot > 0, target_sum / np.maximum(tot, 1e-12), 0.0)
    return np.log1p(X * sf[:, None]).astype(np.float32)
