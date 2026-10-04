"""Exposure mappings (docs/estimands.md, Section 2).

For a recipient cell i and target gene g, the exposure is the vector of counts of
g-perturbed cells in each distance ring around i. Rings are defined by ``bins_um``
(edges, increasing, starting at 0). Everything is computed within samples from the sparse
pairwise-distance matrix produced by :func:`spatialspill.graphs.pairwise_distances_within`.

The ring adjacency matrices depend only on geometry, so they are computed once; the
label-dependent part (:meth:`Exposure.recompute`) is cheap and is what the permutation null
re-evaluates.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import pairwise

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from spatialspill.graphs import pairwise_distances_within

NTC_PREFIX = "NTC:"


def ring_indicator_matrices(D: sp.csr_matrix, bins_um: list[float]) -> list[sp.csr_matrix]:
    """Split a sparse distance matrix into 0/1 ring adjacency matrices (one per bin)."""
    Dc = D.tocoo()
    out = []
    for lo, hi in pairwise(bins_um):
        m = (Dc.data > lo) & (Dc.data <= hi)
        R = sp.coo_matrix((np.ones(int(m.sum())), (Dc.row[m], Dc.col[m])), shape=Dc.shape).tocsr()
        out.append(R)
    return out


def exposure_labels(obs: pd.DataFrame) -> np.ndarray:
    """Per-cell label: target gene, ``"NTC:<guide>"`` for controls, ``"none"`` otherwise."""
    label = obs["target"].astype(str).to_numpy().copy()
    ntc = obs["is_ntc"].to_numpy()
    label[ntc] = np.char.add(NTC_PREFIX, obs["guide"].astype(str).to_numpy()[ntc])
    return label


@dataclass
class Exposure:
    """Ring geometry plus label-dependent ring counts for every (cell, target) pair."""

    targets: list[str]
    bins_um: list[float]
    rings: list[sp.csr_matrix]
    labels: np.ndarray
    is_perturbed_label: dict[str, bool]
    counts: list[sp.csr_matrix] = field(default_factory=list)
    any_perturbed_within: np.ndarray = field(default_factory=lambda: np.empty(0))
    n_neighbors_within: np.ndarray = field(default_factory=lambda: np.empty(0))

    def __post_init__(self) -> None:
        W = self.rings[0].copy()
        for R in self.rings[1:]:
            W = W + R
        self._W = W.tocsr()
        self.n_neighbors_within = np.asarray(self._W.sum(axis=1)).ravel()
        self.recompute(self.labels)

    @property
    def n_bins(self) -> int:
        return len(self.bins_um) - 1

    @property
    def n_cells(self) -> int:
        return self.rings[0].shape[0]

    def label_matrix(self, labels: np.ndarray) -> sp.csr_matrix:
        tidx = {t: k for k, t in enumerate(self.targets)}
        cols = np.array([tidx.get(t, -1) for t in labels])
        keep = cols >= 0
        return sp.coo_matrix(
            (np.ones(int(keep.sum())), (np.flatnonzero(keep), cols[keep])),
            shape=(len(labels), len(self.targets)),
        ).tocsr()

    def recompute(self, labels: np.ndarray) -> Exposure:
        """Recompute ring counts for a new label vector (same geometry)."""
        self.labels = labels
        Z = self.label_matrix(labels)
        self.counts = [(R @ Z).tocsr() for R in self.rings]
        is_pert = np.array([self.is_perturbed_label.get(t, False) for t in labels], dtype=float)
        self.any_perturbed_within = np.asarray(self._W @ is_pert).ravel()
        return self

    def ring_count(self, target: str, b: int) -> np.ndarray:
        k = self.targets.index(target)
        return np.asarray(self.counts[b][:, k].todense()).ravel()

    def total_count_matrix(self) -> sp.csr_matrix:
        tot = self.counts[0].copy()
        for C in self.counts[1:]:
            tot = tot + C
        return tot.tocsr()


def compute_exposure(
    adata: AnnData,
    bins_um: list[float],
    targets: list[str] | None = None,
    sample_key: str = "sample",
    include_ntc: bool = True,
) -> Exposure:
    """Ring counts of each target around each cell.

    NTC guides are included as pseudo-targets (``"NTC:<guide>"``) by default so they can be run
    through every estimator as the mandatory calibration test.
    """
    bins = [float(b) for b in bins_um]
    if bins[0] != 0 or any(np.diff(bins) <= 0):
        raise ValueError("bins_um must start at 0 and be strictly increasing")
    D = pairwise_distances_within(adata, bins[-1], sample_key=sample_key)
    rings = ring_indicator_matrices(D, bins)
    labels = exposure_labels(adata.obs)
    is_pert = adata.obs["is_perturbed"].to_numpy()
    pert_map: dict[str, bool] = {}
    for lab, p in zip(labels, is_pert):
        pert_map[lab] = bool(p) or pert_map.get(lab, False)
    if targets is None:
        targets = sorted(
            t
            for t in np.unique(labels)
            if t != "none" and (include_ntc or not t.startswith(NTC_PREFIX))
        )
    return Exposure(list(targets), bins, rings, labels, pert_map)


def local_density(adata: AnnData, radius_um: float, sample_key: str = "sample") -> np.ndarray:
    """Number of cells within ``radius_um`` (excluding self), per cell."""
    D = pairwise_distances_within(adata, radius_um, sample_key=sample_key)
    return np.asarray((D > 0).sum(axis=1)).ravel().astype(float)
