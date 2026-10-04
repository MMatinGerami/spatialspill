from spatialspill.graphs import GraphConfig, build_graph, degree, pairwise_distances_within


def test_delaunay_symmetric_and_within_sample(toy):
    build_graph(toy, GraphConfig(kind="delaunay", max_edge_um=40))
    A = toy.obsp["spatial_connectivities"]
    assert (A != A.T).nnz == 0
    D = toy.obsp["spatial_distances"]
    assert D.max() <= 40 + 1e-9
    r, c = A.nonzero()
    s = toy.obs["sample"].to_numpy()
    assert (s[r] == s[c]).all()
    assert degree(toy).sum() == A.nnz


def test_radius_matches_pairwise(toy):
    build_graph(toy, GraphConfig(kind="radius", radius_um=20))
    A = toy.obsp["spatial_connectivities"]
    P = pairwise_distances_within(toy, 20)
    assert A.nnz == P.nnz


def test_knn_degree_at_least_k(toy):
    build_graph(toy, GraphConfig(kind="knn", n_neighs=5))
    assert (degree(toy) >= 5).all()


def test_tiny_sample_complete_graph():
    from tests.conftest import make_toy

    a = make_toy(n=3, n_samples=1)
    build_graph(a, GraphConfig(kind="delaunay", max_edge_um=1e9))
    assert a.obsp["spatial_connectivities"].nnz == 6


def test_unknown_kind(toy):
    import pytest

    with pytest.raises(ValueError):
        build_graph(toy, GraphConfig(kind="foo"))
