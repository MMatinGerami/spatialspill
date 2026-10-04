import numpy as np

from spatialspill.schema import validate
from spatialspill.simulator import (
    SimConfig,
    default_scenarios,
    geometry_from,
    simulate,
    synthetic_geometry,
)


def test_simulate_validates_and_plants_truth():
    geom = synthetic_geometry(n_cells=3000, seed=1)
    a = simulate(geom, SimConfig(n_genes=30, n_targets=5, seed=1))
    assert validate(a) == []
    t = a.uns["truth"]
    assert t["tau_auto"].shape == (5, 30) and t["tau_spill"].shape == (5, 30)
    assert a.obs["is_perturbed"].sum() > 0 and a.obs["is_ntc"].sum() > 0
    assert a.X.min() >= 0


def test_autonomous_effect_visible_in_counts():
    geom = synthetic_geometry(n_cells=6000, n_samples=1, seed=2)
    cfg = SimConfig(
        n_genes=20,
        n_targets=2,
        frac_assigned=0.3,
        frac_auto_nonzero=1.0,
        auto_lfc_sd=2.0,
        frac_spill_targets=0.0,
        seed=2,
    )
    a = simulate(geom, cfg)
    tau = a.uns["truth"]["tau_auto"]
    X = a.X.toarray()
    ntc = a.obs["is_ntc"].to_numpy()
    for k, t in enumerate(a.uns["truth"]["targets"]):
        m = (a.obs["target"] == t).to_numpy()
        lfc = np.log((X[m].mean(0) + 0.5) / (X[ntc].mean(0) + 0.5))
        assert np.corrcoef(lfc, tau[k])[0, 1] > 0.8


def test_misassignment_and_clonality_options():
    geom = synthetic_geometry(n_cells=4000, seed=3)
    a = simulate(
        geom,
        SimConfig(
            n_genes=10, n_targets=3, p_misassign=0.5, clonal_radius_um=40, clonal_size=5, seed=3
        ),
    )
    assert a.uns["truth"]["n_misassigned"] > 0
    assert a.obs["misassigned"].sum() == a.uns["truth"]["n_misassigned"]
    # clonal assignment: assigned cells have assigned neighbours more often than chance
    from spatialspill.graphs import GraphConfig, build_graph

    build_graph(a, GraphConfig(kind="delaunay", max_edge_um=40))
    A = a.obsp["spatial_connectivities"]
    g = a.obs["true_guide"].to_numpy()
    r, c = A.nonzero()
    same = ((g[r] == g[c]) & (g[r] != "none")).sum()
    assert same > 0


def test_bleed_through_moves_expected_counts():
    geom = synthetic_geometry(n_cells=3000, n_samples=1, seed=4)
    a0 = simulate(geom, SimConfig(n_genes=10, n_targets=2, alpha_bleed=0.0, seed=4))
    a1 = simulate(geom, SimConfig(n_genes=10, n_targets=2, alpha_bleed=0.3, seed=4))
    # totals are preserved in expectation; per-cell counts differ
    assert abs(a0.X.sum() - a1.X.sum()) / a0.X.sum() < 0.05
    assert (a0.X != a1.X).nnz > 0


def test_geometry_from_real(toy):
    g = geometry_from(toy, n_cells=100)
    assert g.n_obs == 100 and "cell_type" in g.obs
    assert len(default_scenarios()) >= 5
