import numpy as np

from spatialspill.multiple_testing import bh_fdr
from spatialspill.permutation import (
    permutation_null,
    permutation_pvalue,
    permute_within_strata,
    strata_codes,
)


def test_permute_preserves_strata_counts(toy):
    strata = strata_codes(toy.obs)
    labels = toy.obs["target"].to_numpy().astype(object)
    rng = np.random.default_rng(1)
    perm = permute_within_strata(labels, strata, rng)
    for s in np.unique(strata):
        m = strata == s
        assert sorted(labels[m]) == sorted(perm[m])
    assert (perm != labels).any()


def test_no_strata_is_global_shuffle(toy):
    strata = strata_codes(toy.obs, keys=())
    assert (strata == 0).all()


def test_permutation_pvalue_uniform_under_null():
    rng = np.random.default_rng(0)
    labels = rng.integers(0, 2, 500)
    y = rng.normal(size=500)
    strata = np.zeros(500, dtype=int)

    def stat(lab):
        return np.array([y[lab == 1].mean() - y[lab == 0].mean()])

    null = np.stack(list(permutation_null(labels, strata, stat, n_perm=199, seed=0)))
    p = permutation_pvalue(stat(labels), null)
    assert 0 < p[0] <= 1
    assert null.shape == (199, 1)


def test_pvalue_floor():
    null = np.zeros((99, 1))
    p = permutation_pvalue(np.array([10.0]), null)
    assert np.isclose(p[0], 1 / 100)


def test_bh_passthrough_nan():
    p = np.array([0.01, np.nan, 0.04, 0.5])
    q = bh_fdr(p)
    assert np.isnan(q[1])
    assert np.all(q[~np.isnan(q)] >= p[~np.isnan(p)])
