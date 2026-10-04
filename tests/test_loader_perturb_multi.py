import numpy as np
import pytest

from spatialspill.loaders.perturb_multi import recover_counts


def test_recover_counts_inverts_normalisation():
    rng = np.random.default_rng(0)
    counts = rng.poisson(0.7, size=(50, 40)).astype(float)
    counts[:, 0] += 1  # every cell has at least one count of exactly 1
    tot = counts.sum(1)
    ln = np.log1p(counts * 93.0 / tot[:, None])
    assert np.array_equal(recover_counts(ln), counts)


def test_recover_counts_common_factor_is_unidentifiable():
    counts = np.array([[2.0, 4.0, 0.0, 6.0]])
    ln = np.log1p(counts * 93.0 / counts.sum(1)[:, None])
    assert np.array_equal(recover_counts(ln), counts / 2)


def test_recover_counts_rejects_non_integer():
    bad = np.log1p(np.array([[1.0, np.sqrt(2.0), 0.0]]))
    with pytest.raises(ValueError):
        recover_counts(bad)
