"""Stratified permutation null for guide labels.

Guide identity is shuffled within strata (default: sample x cell type), which preserves the
spatial arrangement, cell-type composition and per-stratum guide frequencies while breaking
any link between a cell's own guide and its position. Under the null of no autonomous and no
spillover effect, any statistic computed from permuted labels has the same distribution as
the observed one.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator

import numpy as np
import pandas as pd


def strata_codes(obs: pd.DataFrame, keys: tuple[str, ...] = ("sample", "cell_type")) -> np.ndarray:
    """Integer code per cell identifying its stratum."""
    if not keys:
        return np.zeros(len(obs), dtype=int)
    combined = obs[list(keys)].astype(str).agg("|".join, axis=1)
    return pd.factorize(combined)[0]


def permute_within_strata(
    labels: np.ndarray, strata: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Return a copy of ``labels`` permuted independently within each stratum."""
    out = labels.copy()
    order = np.argsort(strata, kind="stable")
    sorted_strata = strata[order]
    bounds = np.flatnonzero(np.diff(sorted_strata)) + 1
    for block in np.split(order, bounds):
        if len(block) > 1:
            out[block] = labels[rng.permutation(block)]
    return out


def permutation_null(
    labels: np.ndarray,
    strata: np.ndarray,
    statistic: Callable[[np.ndarray], np.ndarray],
    n_perm: int = 200,
    seed: int = 0,
) -> Iterator[np.ndarray]:
    """Yield ``statistic(permuted_labels)`` for ``n_perm`` stratified permutations."""
    rng = np.random.default_rng(seed)
    for _ in range(n_perm):
        yield statistic(permute_within_strata(labels, strata, rng))


def permutation_pvalue(
    observed: np.ndarray, null: np.ndarray, two_sided: bool = True
) -> np.ndarray:
    """Permutation p-value with the +1 correction (Phipson and Smyth 2010).

    ``null`` has shape (n_perm, ...) matching ``observed``.
    """
    observed = np.asarray(observed)
    null = np.asarray(null)
    n = null.shape[0]
    if two_sided:
        exceed = (np.abs(null) >= np.abs(observed)[None, ...]).sum(axis=0)
    else:
        exceed = (null >= observed[None, ...]).sum(axis=0)
    return (exceed + 1) / (n + 1)
