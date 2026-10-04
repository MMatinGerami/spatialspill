import numpy as np

from spatialspill.estimators import E3DoublyRobust
from spatialspill.estimators.e3_dr import exposure_probabilities
from spatialspill.estimators.groups import GroupConfig
from spatialspill.exposure import compute_exposure
from spatialspill.outcomes import lognorm
from spatialspill.permutation import strata_codes
from tests.conftest import make_toy


def test_exposure_probabilities_are_probabilities(toy):
    exp = compute_exposure(toy, [0, 20, 40])
    strata = strata_codes(toy.obs)
    pa, pr, pc = exposure_probabilities(
        exp,
        toy.obs.is_ntc.to_numpy(),
        toy.obs.is_perturbed.to_numpy(),
        strata,
        GroupConfig(control_policy="unperturbed", clean_controls=False),
        n_draws=20,
        seed=0,
    )
    for P in [pa, *pr, pc]:
        assert P.min() >= 0 and P.max() <= 1
    assert pc.mean() > 0


def test_e3_recovers_autonomous_effect():
    a = make_toy(n=1500, n_genes=4, n_samples=2, seed=8, field_um=1500.0)
    m = (a.obs["target"] == "GENE_A").to_numpy()
    a.X[m, 0] += 20
    exp = compute_exposure(a, [0, 15, 30])
    Y = lognorm(a)
    est = E3DoublyRobust(
        n_draws=30,
        n_boot=50,
        min_cells=3,
        groups=GroupConfig(control_policy="unperturbed", clean_controls=False),
    )
    df = est.fit(a, exp, Y, list(a.var_names)).df
    row = df[
        (df.target == "GENE_A")
        & (df.kind == "autonomous")
        & (df.cell_type == "all")
        & (df.outcome == "G0")
    ].iloc[0]
    assert row["estimate"] > 0.5 and row["pvalue"] < 0.01
    assert np.isfinite(df.loc[df.identified, "se"]).all()
    assert set(df["kind"]) == {"autonomous", "spillover"}
