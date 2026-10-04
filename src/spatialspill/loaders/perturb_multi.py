"""Loader for Perturb-Multi (Saunders et al., Cell 2025; Hugging Face xingjiepan/PerturbMulti).

Source file ``RNA_scaled_crispr_screen_20240615.h5ad`` (14.2 GB): 2,206,191 cells x 209 MERFISH
genes. ``raw.X`` holds log1p(normalised) values, not counts; :func:`recover_counts` inverts
the normalisation exactly (see its docstring). Guide assignment: ``singlet_name`` / ``singlet_gene`` are set
for the 79,179 cells with exactly one barcode above the stringent threshold (``n_thresh3 == 1``);
``control`` is the non-targeting class (4,161 cells over many control guides), which makes
this dataset the best calibration ground in the project. Coordinates ``x``, ``y`` are in
micrometres within a ``batch`` (one liver section each, 11 sections); ``fov`` is the imaging
field of view; ``area`` is the segmentation area in square micrometres.

The file is read with h5py in row chunks so that only the requested batches are held in memory.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from spatialspill.schema import SpatialScreenSchema, edge_distance_from_coords

CONTROL_LABEL = "control"


def _cat(obs: h5py.Group, key: str) -> np.ndarray:
    g = obs[key]
    cats = np.array(
        [c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]], dtype=object
    )
    codes = g["codes"][:]
    out = cats[np.clip(codes, 0, None)]
    out[codes < 0] = ""
    return out


def read_obs(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        df = pd.DataFrame(
            {
                "id": [i.decode() if isinstance(i, bytes) else str(i) for i in obs["id"][:]],
                "fov": obs["fov"][:],
                "x": obs["x"][:],
                "y": obs["y"][:],
                "area": obs["area"][:],
                "batch": _cat(obs, "batch"),
                "cell_type": _cat(obs, "cell_type"),
                "cluster_type": _cat(obs, "cluster_type"),
                "singlet_name": _cat(obs, "singlet_name"),
                "singlet_gene": _cat(obs, "singlet_gene"),
                "n_thresh1": obs["n_thresh1"][:],
                "n_thresh3": obs["n_thresh3"][:],
            }
        )
    return df


def recover_counts(lognorm_block: np.ndarray, tol: float = 1e-6) -> np.ndarray:
    """Invert scanpy's normalize_total + log1p exactly.

    The deposited ``raw.X`` is log1p(count * target_sum / total). For every cell the smallest
    non-zero value corresponds to a count of 1, which gives the cell's size factor; dividing
    by it returns integers (verified to 1e-13 on the deposited data, target_sum = 93).
    Raises if any recovered value is further than ``tol`` from an integer.
    """
    E = np.expm1(lognorm_block.astype(np.float64))
    pos = np.where(E > 0, E, np.inf)
    sf = pos.min(axis=1)
    sf = np.where(np.isfinite(sf), sf, 1.0)
    C = E / sf[:, None]
    R = np.round(C)
    if np.abs(C - R).max() > tol:
        raise ValueError("recovered counts are not integers; the normalisation assumption failed")
    return R.astype(np.float32)


def _read_rows(ds: h5py.Dataset, rows: np.ndarray, chunk: int = 100_000) -> sp.csr_matrix:
    """Read selected rows of a dense HDF5 matrix as sparse float32 counts, streaming in chunks."""
    rows = np.sort(rows)
    blocks: list[sp.csr_matrix] = []
    n = ds.shape[0]
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        sel = rows[(rows >= start) & (rows < stop)]
        if len(sel) == 0:
            continue
        block = ds[start:stop]
        blocks.append(sp.csr_matrix(recover_counts(block[sel - start])))
        del block
    return (
        sp.vstack(blocks).tocsr() if blocks else sp.csr_matrix((0, ds.shape[1]), dtype=np.float32)
    )


def load_perturb_multi(
    raw_dir: str | Path = "data/raw/perturb_multi",
    batches: list[str] | None = None,
    max_cells_per_batch: int | None = None,
    seed: int = 0,
) -> AnnData:
    """Build the schema AnnData for the chosen liver sections (``batches``, strings "0".."10").

    ``max_cells_per_batch`` subsamples unassigned cells uniformly at random (assigned cells are
    always kept) and records the downscaling in ``uns``.
    """
    path = Path(raw_dir) / "RNA_scaled_crispr_screen_20240615.h5ad"
    obs = read_obs(path)
    if batches is not None:
        obs = obs[obs["batch"].isin([str(b) for b in batches])]
    rng = np.random.default_rng(seed)
    if max_cells_per_batch is not None:
        keep_idx: list[np.ndarray] = []
        for _, sub in obs.groupby("batch", sort=False):
            assigned = sub.index[sub["singlet_gene"] != ""].to_numpy()
            other = sub.index[sub["singlet_gene"] == ""].to_numpy()
            n_other = max(0, max_cells_per_batch - len(assigned))
            if len(other) > n_other:
                other = rng.choice(other, size=n_other, replace=False)
            keep_idx.append(np.concatenate([assigned, other]))
        obs = obs.loc[np.sort(np.concatenate(keep_idx))]
    rows = obs.index.to_numpy()
    with h5py.File(path, "r") as f:
        X = _read_rows(f["raw"]["X"], rows)
        var_names = [
            v.decode() if isinstance(v, bytes) else str(v) for v in f["raw"]["var"]["_index"][:]
        ]
    obs = obs.reset_index(drop=True)
    guide = obs["singlet_name"].replace("", "none").astype(str)
    target = obs["singlet_gene"].replace("", "none").astype(str)
    is_ntc = (target == CONTROL_LABEL).to_numpy()
    target = target.where(~is_ntc, "NTC")
    out = pd.DataFrame(
        {
            "guide": guide.to_numpy(),
            "guide_confidence": np.nan,
            "target": target.to_numpy(),
            "is_ntc": is_ntc,
            "is_perturbed": (guide != "none").to_numpy() & ~is_ntc,
            "area": obs["area"].to_numpy(),
            "cell_type": obs["cell_type"].replace("", "unknown").to_numpy(),
            "sample": obs["batch"].astype(str).to_numpy(),
            "batch": obs["batch"].astype(str).to_numpy(),
            "fov": obs["fov"].to_numpy(),
            "cluster_type": obs["cluster_type"].to_numpy(),
            "n_barcodes_lenient": obs["n_thresh1"].to_numpy(),
            "n_barcodes_stringent": obs["n_thresh3"].to_numpy(),
        },
        index=obs["id"].to_numpy(),
    )
    xy = obs[["x", "y"]].to_numpy(dtype=float)
    out["edge_distance"] = edge_distance_from_coords(xy, out["sample"].to_numpy())
    out["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    adata = AnnData(X=X, obs=out)
    adata.var_names = var_names
    adata.obsm["spatial"] = xy
    adata.uns["spatialspill"] = SpatialScreenSchema(
        technology="Perturb-Multi (MERFISH)",
        dataset="HF xingjiepan/PerturbMulti crispr_screen_20240615",
        organism="Mus musculus",
        units="um",
        notes={
            "batches": ",".join(sorted(out["sample"].unique())),
            "max_cells_per_batch": str(max_cells_per_batch),
            "guide_call": "singlet_name (n_thresh3 == 1); control guides are NTC",
            "counts": "recovered exactly from log1p-normalised raw.X (target_sum 93)",
        },
    ).to_uns()
    return adata
