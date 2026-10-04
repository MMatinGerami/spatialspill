from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from spatialspill.loaders.spatial_perturbseq import (
    CHIP_FILES,
    GUIDE_TARGETS,
    assign_guides,
    guide_to_target,
    load_spatial_perturbseq,
)
from spatialspill.schema import NTC_LABEL, validate

RAW = Path("data/raw/spatial_perturbseq")


def _gef(chip: str) -> Path | None:
    for p in (RAW / "extracted" / CHIP_FILES[chip], RAW / CHIP_FILES[chip]):
        if p.exists():
            return p
    return None


needs_raw = pytest.mark.skipif(
    _gef("A03599E2") is None or _gef("B03018A2") is None,
    reason="Spatial Perturb-seq GEF files not downloaded",
)


def test_guide_to_target_mapping() -> None:
    assert len(GUIDE_TARGETS) == 18
    assert guide_to_target("sgrna_dpp5") == "Dpp6"
    assert guide_to_target("sgrna_oligo2") == "Olig2"
    assert guide_to_target("sgrna_msafe") == NTC_LABEL
    with pytest.raises(KeyError):
        guide_to_target("sgrna_unknown")


def test_assign_guides_rules() -> None:
    names = ["sgrna_clu", "sgrna_msafe", "sgrna_trem2"]
    G = np.array([[0, 0, 0], [3, 0, 1], [0, 2, 0], [1, 0, 0], [2, 2, 0]], dtype=np.int32)
    df = assign_guides(G, names, min_guide_umi=1)
    assert list(df["guide"]) == ["none", "multi", "sgrna_msafe", "sgrna_clu", "multi"]
    assert list(df["target"]) == ["none", "none", NTC_LABEL, "Clu", "none"]
    assert list(df["is_ntc"]) == [False, False, True, False, False]
    assert list(df["is_perturbed"]) == [False, False, False, True, False]
    assert np.isnan(df["guide_confidence"].iloc[0])
    assert df["guide_confidence"].iloc[1] == pytest.approx(0.75)
    assert list(df["guide_umis_total"]) == [0, 4, 2, 1, 4]
    # threshold 2: a single-UMI cell becomes "none", the 3+1 cell becomes a clean clu call
    df2 = assign_guides(G, names, min_guide_umi=2)
    assert list(df2["guide"]) == ["none", "sgrna_clu", "sgrna_msafe", "none", "multi"]


@needs_raw
def test_load_small_chip_schema() -> None:
    adata = load_spatial_perturbseq(RAW, chips=["A03599E2"])
    assert validate(adata) == []
    assert adata.n_obs == 65_992
    assert not any(v.startswith("sgrna_") for v in adata.var_names)
    assert "" not in set(adata.var_names)
    gc = adata.obsm["guide_counts"]
    assert gc.shape == (adata.n_obs, 18)
    assert set(gc.columns) == {f"sgrna_{k}" for k in GUIDE_TARGETS}
    assert adata.uns["spatialspill"]["notes"]["guide_features"].count("sgrna_") == 18
    targets = set(adata.obs.loc[adata.obs["is_perturbed"], "target"])
    assert len(targets) >= 10
    assert adata.obs["is_ntc"].sum() > 0
    xy = adata.obsm["spatial"]
    extent = xy.max(axis=0) - xy.min(axis=0)
    assert (extent > 5_000).all() and (extent < 12_000).all()
    assert (adata.obs["area"] > 0).all()
    assert (
        adata.obs["n_umi"].to_numpy()
        == np.asarray(adata.X.sum(axis=1)).ravel() + gc.sum(axis=1).to_numpy()
    ).all()


@needs_raw
def test_load_two_chips_concat() -> None:
    adata = load_spatial_perturbseq(RAW, chips=["A03599E2", "B03018A2"])
    assert validate(adata) == []
    assert adata.n_obs == 65_992 + 123_921
    assert set(adata.obs["sample"]) == {"A03599E2", "B03018A2"}
    assert adata.obsm["guide_counts"].shape[1] == 18
    assert adata.X.min() >= 0
