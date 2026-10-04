"""Name -> loader dispatch used by scripts and configs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
from anndata import AnnData

from spatialspill.loaders.perturb_dbit import load_perturb_dbit
from spatialspill.loaders.perturb_fish import load_perturb_fish_tumor
from spatialspill.loaders.perturb_map import load_perturb_map
from spatialspill.loaders.perturb_multi import load_perturb_multi
from spatialspill.loaders.spatial_perturbseq import load_spatial_perturbseq

BUNDLED = Path(__file__).resolve().parents[3] / "data" / "bundled"


def load_dataset(name: str, **kwargs: Any) -> AnnData:
    if name == "bundled":
        return ad.read_h5ad(BUNDLED / "perturb_dbit_mTSG.LV.2_top300.h5ad")
    if name == "perturb_dbit":
        return load_perturb_dbit(**kwargs)
    if name == "perturb_multi":
        return load_perturb_multi(**kwargs)
    if name == "perturb_map":
        return load_perturb_map(**kwargs)
    if name == "perturb_fish":
        return load_perturb_fish_tumor(**kwargs)
    if name == "spatial_perturbseq":
        return load_spatial_perturbseq(**kwargs)
    raise KeyError(f"unknown dataset {name!r}")
