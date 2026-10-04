from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from spatialspill.loaders.perturb_map import load_perturb_map, phenotype_to_target
from spatialspill.schema import validate

RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "perturb_map"
needs_data = pytest.mark.skipif(
    not (RAW / "extracted").is_dir(), reason="Perturb-map raw data not downloaded"
)


def test_phenotype_to_target():
    assert phenotype_to_target("Tgfbr2_1") == "Tgfbr2"
    assert phenotype_to_target("Tgfbr2_4-1") == "Tgfbr2"
    assert phenotype_to_target("Jak2_1") == "Jak2"
    assert phenotype_to_target("KP_1-1") == "KP_unlabeled"
    assert phenotype_to_target("periphery") == "none"
    assert phenotype_to_target("none") == "none"


@pytest.fixture(scope="module")
def adata():
    return load_perturb_map(RAW)


@needs_data
def test_schema_and_samples(adata):
    assert validate(adata) == []
    assert sorted(adata.obs["sample"].unique()) == ["KP_1", "KP_2", "KP_3", "KP_4"]
    assert adata.obs["batch"].nunique() == 4
    assert adata.var_names.is_unique
    assert adata.X.dtype == np.float32


@needs_data
def test_targets(adata):
    targets = set(adata.obs["target"])
    assert {"Tgfbr2", "Jak2"} <= targets
    assert "NTC" not in targets
    assert not adata.obs["is_ntc"].any()
    kp = adata.obs["target"] == "KP_unlabeled"
    assert kp.any() and not adata.obs.loc[kp, "is_perturbed"].any()
    named = ~adata.obs["target"].isin(["none", "KP_unlabeled"])
    assert (adata.obs["is_perturbed"] == named).all()


@needs_data
def test_spot_geometry_in_um(adata):
    from scipy.spatial import cKDTree

    for s in ["KP_1", "KP_4"]:
        xy = adata.obsm["spatial"][(adata.obs["sample"] == s).to_numpy()]
        d, _ = cKDTree(xy).query(xy, k=2)
        assert abs(np.median(d[:, 1]) - 100.0) < 2.0
