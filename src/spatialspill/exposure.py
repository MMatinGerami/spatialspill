"""Exposure mappings (docs/estimands.md, Section 2).

For a recipient cell i and target gene g, the exposure is the vector of counts of
g-perturbed cells in each distance ring around i. Rings are defined by ``bins_um``
(edges, increasing, starting at 0). Everything is computed within samples from the sparse
pairwise-distance matrix produced by :func:`spatialspill.graphs.pairwise_distances_within`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from spatialspill.graphs import pairwise_distances_within


@dataclass
class Exposure:
    """Ring counts for every (cell, target) pair.

    ``counts[b]`` is an (n_cells x n_targets) sparse matrix of the number of cells with
    target ``targets[k]`` in ring ``b`` around each cell. ``any_perturbed_within`` is the
    number of perturbed neighbours (any targeting guide) within D_max.
    """

    targets: list[str]
    bins_um: list[float]
    counts: list[sp.csr_matrix]
    any_perturbed_within: np.ndarray
    n_neighbors_within: np.ndarray

    @property
    def n_bins(self) -> int:
        return len(self.bins_um) - 1

    def ring_count(self, target: str, b: int) -> np.ndarray:
        k = self.targets.index(target)
        return np.asarray(self.counts[b][:, k].todense()).ravel()

    def total_count(self, target: str) -> np.ndarray:
        return sum(self.ring_count(target, b) for b in range(self.n_bins))


def ring_indicator_matrices(D: sp.csr_matrix, bins_um: list[float]) -> list[sp.csr_matrix]:
    """Split a sparse distance matrix into 0/1 ring adjacency matrices (one per bin)."""
    D = D.tocoo()
    out = []
    for lo, hi in zip(bins_um[:-1], bins_um[1:]):
        m = (D.data > lo) & (D.data <= hi)
        R = sp.coo_matrix((np.ones(m.sum()), (D.row[m], D.col[m])), shape=D.shape).tocsr()
        out.append(R)
    return out


def compute_exposure(
    adata: AnnData,
    bins_um: list[float],
    targets: list[str] | None = None,
    target_key: str = "target",
    sample_key: str = "sample",
    include_ntc: bool = True,
) -> Exposure:
    """Ring counts of each target around each cell.

    NTC guides are included as pseudo-targets by default so they can be run through every
    estimator as the mandatory calibration test. Individual NTC guides are distinguished by
    the ``guide`` column, so a target list entry may be ``"NTC:<guide>"``.
    """
    bins = [float(b) for b in bins_um]
    if bins[0] != 0 or any(np.diff(bins) <= 0):
        raise ValueError("bins_um must start at 0 and be strictly increasing")
    D = pairwise_distances_within(adata, bins[-1], sample_key=sample_key)
    rings = ring_indicator_matrices(D, bins)

    obs = adata.obs
    label = obs[target_key].astype(str).to_numpy().copy()
    if include_ntc:
        ntc = obs["is_ntc"].to_numpy()
        label[ntc] = np.char.add("NTC:", obs["guide"].astype(str).to_numpy()[ntc])
    if targets is None:
        targets = sorted(t for t in np.unique(label) if t not in ("none", "NTC"))
    tidx = {t: k for k, t in enumerate(targets)}
    cols = np.array([tidx.get(t, -1) for t in label])
    keep = cols >= 0
    Z = sp.coo_matrix(
        (np.ones(keep.sum()), (np.flatnonzero(keep), cols[keep])),
        shape=(adata.n_obs, len(targets)),
    ).tocsr()
    counts = [(R @ Z).tocsr() for R in rings]
    is_pert = obs["is_perturbed"].to_numpy().astype(float)
    W = sum(rings)
    any_pert = np.asarray(W @ is_pert).ravel()
    n_within = np.asarray(W.sum(axis=1)).ravel()
    return Exposure(list(targets), bins, counts, any_pert, n_within)


def local_density(adata: AnnData, radius_um: float, sample_key: str = "sample") -> np.ndarray:
    """Number of cells within ``radius_um`` (excluding self), per cell."""
    D = pairwise_distances_within(adata, radius_um, sample_key=sample_key)
    return np.asarray((D > 0).sum(axis=1)).ravel().astype(float)


def exposure_table(exp: Exposure, target: str) -> pd.DataFrame:
    """Long table (cell, ring, count) for one target, used in reports."""
    rows = []
    for b in range(exp.n_bins):
        rows.append(
            pd.DataFrame(
                {
                    "cell": np.arange(exp.counts[b].shape[0]),
                    "ring": b,
                    "count": exp.ring_count(target, b),
                }
            )
        )
    return pd.concat(rows, ignore_index=True)
