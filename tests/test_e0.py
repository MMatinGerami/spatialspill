from spatialspill.estimators import NeighbourTTest, PseudobulkDE
from spatialspill.estimators.groups import GroupConfig
from spatialspill.exposure import compute_exposure
from spatialspill.outcomes import lognorm
from tests.conftest import make_toy


def test_baselines_run_and_recover_autonomous():
    a = make_toy(n=1500, n_genes=5, n_samples=3, seed=9, field_um=1500.0)
    m = (a.obs["target"] == "GENE_A").to_numpy()
    a.X[m, 0] += 20
    exp = compute_exposure(a, [0, 15, 30])
    Y = lognorm(a)
    g = GroupConfig(control_policy="unperturbed", clean_controls=False)
    for est in (NeighbourTTest(min_cells=3, groups=g), PseudobulkDE(min_cells=3, groups=g)):
        df = est.fit(a, exp, Y, list(a.var_names)).df
        row = df[(df.target == "GENE_A") & (df.kind == "autonomous") & (df.outcome == "G0")].iloc[0]
        assert row["estimate"] > 0.5 and row["pvalue"] < 0.01
    dfn = NeighbourTTest(min_cells=3, groups=g).fit(a, exp, Y, list(a.var_names)).df
    assert (dfn.kind == "spillover").any()
    dfp = PseudobulkDE(min_cells=3, groups=g).fit(a, exp, Y, list(a.var_names)).df
    assert not (dfp.kind == "spillover").any()
