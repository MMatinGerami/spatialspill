import numpy as np

from spatialspill.estimators import E4GNN
from spatialspill.estimators.base import COLUMNS
from spatialspill.estimators.groups import GroupConfig
from spatialspill.exposure import compute_exposure
from spatialspill.outcomes import lognorm
from spatialspill.simulator import SimConfig, simulate, synthetic_geometry
from tests.conftest import make_toy


def _toy_fit():
    a = make_toy(n=1500, n_genes=5, n_samples=2, seed=5, field_um=1500.0)
    m = (a.obs["target"] == "GENE_A").to_numpy()
    a.X[m, 0] += 20
    exp = compute_exposure(a, [0, 15, 30])
    Y = lognorm(a)
    est = E4GNN(
        min_cells=3,
        device="cpu",
        groups=GroupConfig(control_policy="unperturbed", clean_controls=False),
    )
    return est.fit(a, exp, Y, list(a.var_names))


def test_e4_recovers_planted_autonomous_effect():
    df = _toy_fit().df
    row = df[
        (df.target == "GENE_A")
        & (df.kind == "autonomous")
        & (df.cell_type == "all")
        & (df.outcome == "G0")
    ].iloc[0]
    assert row["estimate"] > 0.3
    assert row["ci_low"] < row["estimate"] < row["ci_high"]
    assert row["identified"] and row["pvalue"] < 1e-3
    # the target without a planted effect stays small (GENE_A's other genes shift through
    # the size factor, so they are not a null)
    other = df[(df.kind == "autonomous") & (df.cell_type == "all") & (df.target == "GENE_B")]
    assert (other["estimate"].abs() < 0.3).all()


def test_e4_output_schema():
    tab = _toy_fit()
    df = tab.df
    assert list(df.columns) == COLUMNS
    assert set(df["kind"]) == {"autonomous", "spillover"}
    assert set(df.loc[df.kind == "autonomous", "ring"]) == {-1}
    assert set(df.loc[df.kind == "spillover", "ring"]) == {0, 1}
    assert (df["estimator"] == "E4_gnn").all()
    assert {"A", "B", "all"} <= set(df["cell_type"])
    assert np.isfinite(df.loc[df.identified, ["estimate", "se", "pvalue"]].to_numpy()).all()
    q = tab.with_fdr()
    assert q.loc[q["identified"], "qvalue"].notna().all()


def test_e4_recovers_planted_spillover_direction():
    geom = synthetic_geometry(n_cells=4000, n_samples=1, seed=0)
    sim = simulate(
        geom,
        SimConfig(
            n_genes=10,
            n_targets=2,
            frac_assigned=0.3,
            frac_spill_targets=1.0,
            frac_spill_genes=1.0,
            spill_lfc_sd=2.0,
            frac_auto_nonzero=0.0,
            seed=0,
        ),
    )
    exp = compute_exposure(sim, [0, 30, 60])
    Y = lognorm(sim)
    df = (
        E4GNN(
            device="cpu",
            groups=GroupConfig(control_policy="any_other_guide", clean_controls=False),
        )
        .fit(sim, exp, Y, list(sim.var_names))
        .df
    )
    truth = sim.uns["truth"]
    est, tru = [], []
    for k, t in enumerate(truth["targets"]):
        sub = df[
            (df.target == t) & (df.kind == "spillover") & (df.ring == 0) & (df.cell_type == "all")
        ]
        sub = sub.set_index("outcome").reindex(list(sim.var_names))
        est.append(sub["estimate"].to_numpy())
        tru.append(truth["tau_spill"][k])
    est_a, tru_a = np.concatenate(est), np.concatenate(tru)
    ok = np.isfinite(est_a)
    assert ok.sum() >= 10
    r = np.corrcoef(est_a[ok], tru_a[ok])[0, 1]
    print(f"ring-0 spillover vs truth: r = {r:.3f} over {int(ok.sum())} (target, gene) pairs")
    assert r > 0.3, r
