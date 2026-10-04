"""Spatial neighbour graphs: Delaunay and radius graphs, per sample, via squidpy.

Default (justified in docs/estimands.md and reports/audit): Delaunay graph with edges longer
than ``max_edge_um`` pruned, built separately within each sample so no edge crosses fields of
view. A radius graph is the sensitivity alternative. Both return a symmetric sparse adjacency
in ``adata.obsp[key]`` plus a matching distance matrix.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from anndata import AnnData
from scipy.spatial import Delaunay, cKDTree


@dataclass(frozen=True)
class GraphConfig:
    kind: str = "delaunay"  # "delaunay" | "radius" | "knn"
    max_edge_um: float = 50.0  # prune Delaunay edges longer than this
    radius_um: float = 30.0  # for kind == "radius"
    n_neighs: int = 10  # for kind == "knn"
    sample_key: str = "sample"
    key_added: str = "spatial"


def _delaunay_edges(xy: np.ndarray) -> np.ndarray:
    if len(xy) < 4:
        # fall back to complete graph on tiny sets
        idx = np.arange(len(xy))
        i, j = np.meshgrid(idx, idx, indexing="ij")
        m = i < j
        return np.stack([i[m], j[m]], axis=1)
    tri = Delaunay(xy)
    s = tri.simplices
    e = np.concatenate([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]], axis=0)
    e = np.sort(e, axis=1)
    return np.unique(e, axis=0)


def build_graph(adata: AnnData, cfg: GraphConfig | None = None) -> AnnData:
    """Build a per-sample neighbour graph and store it in ``adata.obsp``.

    Writes ``obsp[f"{key}_connectivities"]`` (0/1, symmetric) and
    ``obsp[f"{key}_distances"]`` (micrometres), and ``uns[f"{key}_neighbors"]`` with the
    config, mirroring squidpy's layout so ``squidpy.gr.*`` functions accept it.
    """
    if cfg is None:
        cfg = GraphConfig()
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    n = len(xy)
    samples = adata.obs[cfg.sample_key].to_numpy()
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    dists: list[np.ndarray] = []
    for s in np.unique(samples):
        idx = np.flatnonzero(samples == s)
        pts = xy[idx]
        if len(idx) < 2:
            continue
        if cfg.kind == "delaunay":
            e = _delaunay_edges(pts)
            d = np.linalg.norm(pts[e[:, 0]] - pts[e[:, 1]], axis=1)
            keep = d <= cfg.max_edge_um
            e, d = e[keep], d[keep]
        elif cfg.kind == "radius":
            tree = cKDTree(pts)
            pairs = tree.query_pairs(cfg.radius_um, output_type="ndarray")
            e = pairs
            d = np.linalg.norm(pts[e[:, 0]] - pts[e[:, 1]], axis=1) if len(e) else np.empty(0)
        elif cfg.kind == "knn":
            tree = cKDTree(pts)
            k = min(cfg.n_neighs + 1, len(pts))
            dd, ii = tree.query(pts, k=k)
            src = np.repeat(np.arange(len(pts)), k - 1)
            dst = ii[:, 1:].ravel()
            d = dd[:, 1:].ravel()
            e = np.sort(np.stack([src, dst], axis=1), axis=1)
            e, uidx = np.unique(e, axis=0, return_index=True)
            d = d[uidx]
        else:
            raise ValueError(f"unknown graph kind {cfg.kind!r}")
        rows.append(idx[e[:, 0]])
        cols.append(idx[e[:, 1]])
        dists.append(d)
    if rows:
        r = np.concatenate(rows)
        c = np.concatenate(cols)
        d = np.concatenate(dists)
    else:
        r = c = np.empty(0, dtype=int)
        d = np.empty(0)
    A = sp.coo_matrix((np.ones(len(r)), (r, c)), shape=(n, n))
    D = sp.coo_matrix((d, (r, c)), shape=(n, n))
    A = (A + A.T).tocsr()
    D = (D + D.T).tocsr()
    A.data[:] = 1.0
    key = cfg.key_added
    adata.obsp[f"{key}_connectivities"] = A
    adata.obsp[f"{key}_distances"] = D
    adata.uns[f"{key}_neighbors"] = {
        "connectivities_key": f"{key}_connectivities",
        "distances_key": f"{key}_distances",
        "params": {
            "kind": cfg.kind,
            "max_edge_um": cfg.max_edge_um,
            "radius_um": cfg.radius_um,
            "n_neighs": cfg.n_neighs,
            "sample_key": cfg.sample_key,
        },
    }
    return adata


def pairwise_distances_within(
    adata: AnnData, max_um: float, sample_key: str = "sample"
) -> sp.csr_matrix:
    """Sparse matrix of all pairwise distances <= ``max_um`` within each sample.

    Used for distance-binned estimands; independent of the neighbour-graph definition.
    """
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    n = len(xy)
    samples = adata.obs[sample_key].to_numpy()
    rows, cols, dists = [], [], []
    for s in np.unique(samples):
        idx = np.flatnonzero(samples == s)
        if len(idx) < 2:
            continue
        tree = cKDTree(xy[idx])
        pairs = tree.query_pairs(max_um, output_type="ndarray")
        if len(pairs) == 0:
            continue
        d = np.linalg.norm(xy[idx][pairs[:, 0]] - xy[idx][pairs[:, 1]], axis=1)
        rows.append(idx[pairs[:, 0]])
        cols.append(idx[pairs[:, 1]])
        dists.append(d)
    if not rows:
        return sp.csr_matrix((n, n))
    r = np.concatenate(rows)
    c = np.concatenate(cols)
    d = np.concatenate(dists)
    D = sp.coo_matrix((d, (r, c)), shape=(n, n))
    return (D + D.T).tocsr()


def degree(adata: AnnData, key: str = "spatial") -> np.ndarray:
    A = adata.obsp[f"{key}_connectivities"]
    return np.asarray(A.sum(axis=1)).ravel()
