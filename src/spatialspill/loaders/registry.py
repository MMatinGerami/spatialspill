"""Name -> loader dispatch used by scripts and configs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
from anndata import AnnData

from spatialspill.loaders.perturb_dbit import load_perturb_dbit
from spatialspill.loaders.perturb_multi import load_perturb_multi

BUNDLED = Path(__file__).resolve().parents[3] / "data" / "bundled"


def load_dataset(name: str, **kwargs: Any) -> AnnData:
    if name == "bundled":
        return ad.read_h5ad(BUNDLED / "perturb_dbit_mTSG.LV.2_top300.h5ad")
    if name == "perturb_dbit":
        return load_perturb_dbit(**kwargs)
    if name == "perturb_multi":
        return load_perturb_multi(**kwargs)
    raise KeyError(f"unknown dataset {name!r}")
