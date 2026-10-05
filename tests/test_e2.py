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


def test_e2_spatial_basis_runs_and_absorbs_smooth_field():
    a = make_toy(n=2000, n_genes=4, n_samples=2, seed=12, field_um=1500.0)
    # a smooth spatial field on G0 that is unrelated to guides
    xy = a.obsm["spatial"]
    a.X[:, 0] = a.X[:, 0] + np.round(5 * (1 + np.sin(xy[:, 0] / 300) * np.cos(xy[:, 1] / 300)))
    exp = compute_exposure(a, [0, 20, 40])
    Y = lognorm(a)
    g = GroupConfig(control_policy="unperturbed", clean_controls=False)
    plain = E2GLM(min_cells=3, groups=g).fit(a, exp, Y, list(a.var_names)).df
    spatial = E2GLM(min_cells=3, groups=g, spatial_basis=25).fit(a, exp, Y, list(a.var_names)).df
    assert spatial["identified"].sum() > 0
    s0 = spatial[(spatial.outcome == "G0") & spatial.identified]["se"].median()
    p0 = plain[(plain.outcome == "G0") & plain.identified]["se"].median()
    assert s0 <= p0 * 1.05  # absorbing the field should not inflate uncertainty


def test_e2_ring_coefficients_ignore_own_cell_neighbourhoods():
    # two targets whose cells are clustered (clone-like); GENE_A cells have a huge, homogeneous
    # shift in G0; recipients are unaffected, so the ring coefficient must be near zero.
    a = make_toy(n=3000, n_genes=4, n_samples=1, seed=21, field_um=2000.0, p_perturbed=0.3)
    xy = a.obsm["spatial"]
    clone = (xy[:, 0] < 300) & (xy[:, 1] < 300)
    a.obs.loc[clone, "target"] = "GENE_A"
    a.obs.loc[clone, "guide"] = "GENE_A_g1"
    a.obs.loc[clone, "is_ntc"] = False
    a.obs.loc[clone, "is_perturbed"] = True
    a.X[clone, 0] += 50
    exp = compute_exposure(a, [0, 20, 40])
    Y = lognorm(a)
    df = (
        E2GLM(min_cells=3, groups=GroupConfig(control_policy="unperturbed", clean_controls=False))
        .fit(a, exp, Y, list(a.var_names))
        .df
    )
    r = df[
        (df.target == "GENE_A")
        & (df.kind == "spillover")
        & (df.cell_type == "all")
        & (df.outcome == "G0")
        & df.identified
    ]
    assert len(r) > 0
    assert (r["pvalue"] > 1e-4).all()
