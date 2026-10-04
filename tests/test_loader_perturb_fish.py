from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from spatialspill.loaders.perturb_fish import (
    TARGETS,
    assign_guides,
    guide_names,
    load_perturb_fish_tumor,
)
from spatialspill.schema import validate

RAW = Path("data/raw/perturb_fish")
FINAL = RAW / "extras" / "tumors" / "processed" / "finaltables"
REQUIRED = (
    "merfishcounttable.csv",
    "coordinates.csv",
    "allcellsPerturbationTable.csv",
    "tumorMerfish.csv",
)

needs_raw = pytest.mark.skipif(
    not all((FINAL / f).exists() for f in REQUIRED),
    reason="Perturb-FISH raw tables not downloaded",
)


def test_guide_names_layout() -> None:
    names = guide_names()
    assert len(names) == 77
    assert names[0] == "CD14_g1" and names[1] == "CD14_g2"
    assert names[69] == "TRAM1_g2"
    assert names[70:74] == [f"Control_g{k}" for k in range(1, 5)]
    assert len(TARGETS) == 35 and "LBP" in TARGETS


def test_assign_guides_rules() -> None:
    P = np.zeros((5, 77), dtype=np.int8)
    P[0, 0] = 1  # CD14_g1
    P[1, 71] = 1  # Control_g2
    P[2, 0] = 1
    P[2, 2] = 1  # two guides -> multi
    P[3, 6] = 1  # noisy column 7 -> zeroed
    P[4, 74] = 1  # unused column 75 -> zeroed
    df = assign_guides(P)
    assert df["guide"].tolist() == ["CD14_g1", "Control_g2", "multi", "none", "none"]
    assert df["target"].tolist() == ["CD14", "NTC", "none", "none", "none"]
    assert df["is_ntc"].tolist() == [False, True, False, False, False]
    assert df["is_perturbed"].tolist() == [True, False, False, False, False]
    assert df["n_guides"].tolist() == [1, 1, 2, 0, 0]


@needs_raw
def test_load_perturb_fish_tumor_schema() -> None:
    adata = load_perturb_fish_tumor(RAW)
    assert validate(adata) == []
    assert adata.n_obs == 187_215
    assert adata.n_vars == 500
    assert not adata.var["is_blank"].any()
    xy = adata.obsm["spatial"]
    span = xy.max(axis=0) - xy.min(axis=0)
    # A single MERSCOPE section: a few millimetres on each side.
    assert 2_000 < span[0] < 15_000 and 2_000 < span[1] < 15_000
    assert xy.min() >= 0
    obs = adata.obs
    targets = obs.loc[obs["is_perturbed"], "target"].unique()
    assert len(targets) >= 20
    assert obs["is_ntc"].sum() > 1_000
    assert (obs.loc[obs["guide"] == "multi", "n_guides"] >= 2).all()
    assert not obs.loc[obs["guide"] == "multi", "is_perturbed"].any()
    assert (obs["cell_type"] == "T cell").sum() > 10_000
    assert (obs["total_counts"] == np.asarray(adata.X.sum(axis=1)).ravel()).all()
