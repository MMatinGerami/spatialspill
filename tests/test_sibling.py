import numpy as np

from spatialspill.exposure import compute_exposure
from spatialspill.outcomes import lognorm
from spatialspill.sibling import estimate_pi, sibling_corrected_spillover, sibling_table
from spatialspill.simulator import SimConfig, simulate, synthetic_geometry


def test_estimate_pi_recovers_mixture():
    rng = np.random.default_rng(0)
    a = rng.normal(0, 1, 50)
    s = rng.normal(0, 0.3, 50)
    u = 0.4 * a + 0.6 * s
    assert abs(estimate_pi(u, a) - 0.4) < 0.1


def test_sibling_correction_reduces_autonomous_contamination():
    geom = synthetic_geometry(n_cells=8000, n_samples=1, seed=7)
    cfg = SimConfig(
        n_genes=30,
        n_targets=3,
        frac_assigned=0.3,
        frac_auto_nonzero=0.5,
        auto_lfc_sd=1.5,
        frac_spill_targets=0.0,
        p_misassign=0.0,
        clonal_radius_um=40.0,
        clonal_size=6,
        seed=7,
    )
    sim = simulate(geom, cfg)
    # hide 60% of the guide calls so that unassigned neighbours are mostly undetected siblings
    rng = np.random.default_rng(1)
    assigned = sim.obs["guide"] != "none"
    hide = assigned.to_numpy() & (rng.random(sim.n_obs) < 0.6)
    sim.obs.loc[hide, ["guide", "target"]] = "none"
    sim.obs.loc[hide, "is_ntc"] = False
    sim.obs.loc[hide, "is_perturbed"] = False
    exp = compute_exposure(sim, [0, 15, 30, 60])
    Y = lognorm(sim)
    tgt = sim.uns["truth"]["targets"][0]
    est = sibling_corrected_spillover(sim, exp, Y, tgt, reference="ntc", n_boot=50)
    assert est is not None
    a = est.autonomous
    # uncorrected unassigned-neighbour profile correlates with the autonomous profile
    r_unc = np.corrcoef(est.uncorrected, a)[0, 1]
    r_cor = np.corrcoef(est.corrected, a)[0, 1]
    assert est.pi > 0.05  # siblings are a minority of unassigned neighbours at 60 um
    assert abs(r_cor) < abs(r_unc)
    tab = sibling_table(sim, exp, Y, list(sim.var_names), n_boot=20)
    assert {"pi", "corrected", "uncorrected"} <= set(tab.columns)
