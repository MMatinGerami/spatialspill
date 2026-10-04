import numpy as np

from spatialspill.estimators import E2GLM
from spatialspill.estimators.groups import GroupConfig
from spatialspill.exposure import compute_exposure
from spatialspill.outcomes import lognorm
from tests.conftest import make_toy


def test_e2_runs_and_recovers_autonomous_effect():
    a = make_toy(n=1500, n_genes=5, n_samples=2, seed=5, field_um=1500.0)
    m = (a.obs["target"] == "GENE_A").to_numpy()
    a.X[m, 0] += 20
    exp = compute_exposure(a, [0, 15, 30])
    Y = lognorm(a)
    tab = E2GLM(
        min_cells=3, groups=GroupConfig(control_policy="unperturbed", clean_controls=False)
    ).fit(a, exp, Y, list(a.var_names))
    df = tab.df
    assert set(df["kind"]) == {"autonomous", "spillover"}
    row = df[
        (df.target == "GENE_A")
        & (df.kind == "autonomous")
        & (df.cell_type == "all")
        & (df.outcome == "G0")
    ].iloc[0]
    assert row["estimate"] > 0.5 and row["pvalue"] < 1e-3
    assert row["ci_low"] < row["estimate"] < row["ci_high"]
    q = tab.with_fdr()
    assert q.loc[q["identified"], "qvalue"].notna().all()


def test_e2_null_pvalues_roughly_uniform():
    a = make_toy(n=3000, n_genes=20, n_samples=3, seed=6, field_um=2000.0)
    exp = compute_exposure(a, [0, 20, 40])
    Y = lognorm(a)
    df = (
        E2GLM(min_cells=3, groups=GroupConfig(control_policy="unperturbed", clean_controls=False))
        .fit(a, exp, Y, list(a.var_names))
        .df
    )
    p = df.loc[df["identified"], "pvalue"].to_numpy()
    assert 0.01 < (p < 0.05).mean() < 0.12


def test_e2_without_samples_uses_hc1():
    a = make_toy(n=800, n_genes=4, n_samples=1, seed=7, field_um=1200.0)
    exp = compute_exposure(a, [0, 20])
    df = (
        E2GLM(min_cells=3, groups=GroupConfig(control_policy="unperturbed", clean_controls=False))
        .fit(a, exp, lognorm(a), list(a.var_names))
        .df
    )
    assert np.isfinite(df.loc[df["identified"], "se"]).all()
