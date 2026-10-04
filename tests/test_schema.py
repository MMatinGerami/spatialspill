import numpy as np
import pytest

from spatialspill.schema import SchemaError, edge_distance_from_coords, validate


def test_toy_validates(toy):
    assert validate(toy) == []


def test_missing_column_raises(toy):
    bad = toy.copy()
    del bad.obs["area"]
    with pytest.raises(SchemaError, match="missing obs columns"):
        validate(bad)
    assert validate(bad, strict=False)


def test_ntc_and_perturbed_conflict(toy):
    bad = toy.copy()
    bad.obs["is_perturbed"] = True
    probs = validate(bad, strict=False)
    assert any("both NTC and perturbed" in p for p in probs)


def test_negative_counts(toy):
    bad = toy.copy()
    bad.X[0, 0] = -1
    assert any("non-negative" in p for p in validate(bad, strict=False))


def test_edge_distance_bbox():
    xy = np.array([[0.0, 0.0], [10.0, 10.0], [5.0, 5.0]])
    s = np.array(["a", "a", "a"])
    d = edge_distance_from_coords(xy, s)
    assert d[0] == 0 and d[1] == 0 and d[2] == 5
