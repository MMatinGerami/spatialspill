import numpy as np

from spatialspill.loaders.perturb_multi import recover_counts


def test_recover_counts_inverts_normalisation():
    rng = np.random.default_rng(0)
    counts = rng.poisson(0.7, size=(50, 40)).astype(float)
    counts[:, 0] += 1  # every cell has at least one count-1-or-more gene
    tot = counts.sum(1)
    ln = np.log1p(counts * 93.0 / tot[:, None])
    rec = recover_counts(ln)
    assert np.array_equal(rec, counts)


def test_recover_counts_rejects_non_integer():
    import pytest

    bad = np.log1p(np.array([[0.3, 0.5, 0.0]]))
    with pytest.raises(ValueError):
        recover_counts(bad, tol=1e-9)
