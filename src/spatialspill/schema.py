"""Common AnnData/MuData schema for spatial perturbation screens.

Every loader produces an :class:`anndata.AnnData` (one modality, ``rna``) wrapped in a
:class:`mudata.MuData` with these conventions:

``adata.X``            raw counts (cells x genes), sparse or dense, non-negative.
``adata.obsm["spatial"]`` (n, 2) float coordinates in micrometres.
``adata.obs`` columns listed in :data:`OBS_REQUIRED`:

- ``guide``: guide identifier as called by the original study (string; ``"none"`` if no call).
- ``guide_confidence``: float in [0, 1] or NaN when the study provides no confidence.
- ``target``: target gene symbol, ``"NTC"`` for non-targeting controls, ``"none"`` if no call.
- ``is_ntc``: bool, True for non-targeting / safe-harbour controls.
- ``is_perturbed``: bool, True when a targeting guide was called.
- ``area``: segmentation area in square micrometres (NaN if unavailable).
- ``cell_type``: string label (``"unknown"`` allowed).
- ``sample``: biological sample / field of view / slide identifier.
- ``batch``: processing batch (may equal ``sample``).
- ``edge_distance``: distance in micrometres to the tissue or field boundary.

``adata.uns["spatialspill"]`` holds dataset-level metadata: ``technology``, ``dataset``,
``organism``, ``units`` and the loader version.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from anndata import AnnData

OBS_REQUIRED: tuple[str, ...] = (
    "guide",
    "guide_confidence",
    "target",
    "is_ntc",
    "is_perturbed",
    "area",
    "cell_type",
    "sample",
    "batch",
    "edge_distance",
)

UNS_REQUIRED: tuple[str, ...] = ("technology", "dataset", "organism", "units")

NTC_LABEL = "NTC"
NONE_LABEL = "none"


@dataclass
class SpatialScreenSchema:
    """Metadata stored under ``adata.uns["spatialspill"]``."""

    technology: str
    dataset: str
    organism: str
    units: str = "um"
    loader_version: str = "0.1.0"
    notes: dict[str, str] = field(default_factory=dict)

    def to_uns(self) -> dict[str, object]:
        return {
            "technology": self.technology,
            "dataset": self.dataset,
            "organism": self.organism,
            "units": self.units,
            "loader_version": self.loader_version,
            "notes": dict(self.notes),
        }


class SchemaError(ValueError):
    """Raised when an AnnData does not satisfy the spatialspill schema."""


def validate(adata: AnnData, *, strict: bool = True) -> list[str]:
    """Check an AnnData against the schema.

    Returns a list of problems. If ``strict`` is True, raises :class:`SchemaError`
    when the list is non-empty.
    """
    problems: list[str] = []
    missing = [c for c in OBS_REQUIRED if c not in adata.obs.columns]
    if missing:
        problems.append(f"missing obs columns: {missing}")
    if "spatial" not in adata.obsm:
        problems.append("missing obsm['spatial']")
    else:
        xy = np.asarray(adata.obsm["spatial"])
        if xy.ndim != 2 or xy.shape[1] != 2:
            problems.append(f"obsm['spatial'] must be (n, 2), got {xy.shape}")
        elif not np.isfinite(xy).all():
            problems.append("obsm['spatial'] contains non-finite values")
    uns = adata.uns.get("spatialspill")
    if not isinstance(uns, dict):
        problems.append("missing uns['spatialspill']")
    else:
        missing_uns = [k for k in UNS_REQUIRED if k not in uns]
        if missing_uns:
            problems.append(f"missing uns['spatialspill'] keys: {missing_uns}")
    if not missing:
        obs = adata.obs
        if obs["is_ntc"].dtype != bool:
            problems.append("is_ntc must be bool")
        if obs["is_perturbed"].dtype != bool:
            problems.append("is_perturbed must be bool")
        if (obs["is_ntc"] & obs["is_perturbed"]).any():
            problems.append("cells cannot be both NTC and perturbed")
        ntc_targets = obs.loc[obs["is_ntc"], "target"].astype(str)
        if not (ntc_targets == NTC_LABEL).all():
            problems.append(f"NTC cells must have target == {NTC_LABEL!r}")
        conf = pd.to_numeric(obs["guide_confidence"], errors="coerce")
        if ((conf < 0) | (conf > 1)).any():
            problems.append("guide_confidence must lie in [0, 1] or be NaN")
        if (pd.to_numeric(obs["area"], errors="coerce") < 0).any():
            problems.append("area must be non-negative")
        if (pd.to_numeric(obs["edge_distance"], errors="coerce") < 0).any():
            problems.append("edge_distance must be non-negative")
    X = adata.X
    if X is not None:
        mn = X.min() if hasattr(X, "min") else np.min(X)
        if mn < 0:
            problems.append("X must hold non-negative raw counts")
    if strict and problems:
        raise SchemaError("; ".join(problems))
    return problems


def edge_distance_from_coords(xy: np.ndarray, sample: np.ndarray | pd.Series) -> np.ndarray:
    """Distance from each cell to the bounding box of its sample.

    A bounding-box proxy is used when no tissue mask is available. Loaders that have a
    tissue mask or convex hull should override this.
    """
    xy = np.asarray(xy, dtype=float)
    sample = np.asarray(sample)
    out = np.empty(len(xy), dtype=float)
    for s in np.unique(sample):
        m = sample == s
        pts = xy[m]
        lo = pts.min(axis=0)
        hi = pts.max(axis=0)
        d = np.minimum(pts - lo, hi - pts).min(axis=1)
        out[m] = d
    return out
