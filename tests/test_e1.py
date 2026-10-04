import numpy as np

from spatialspill.estimators import E1Stratified
from spatialspill.estimators.groups import GroupConfig, group_masks
from spatialspill.exposure import compute_exposure
from spatialspill.outcomes import lognorm
from tests.conftest import make_toy


def test_group_masks_are_disjoint_and_consistent(toy):
    exp = compute_exposure(toy, [0, 20, 40])
    T, rings, C = group_masks(
        exp, toy.obs["is_ntc"].to_numpy(), toy.obs["is_perturbed"].to_numpy(), GroupConfig()
    )
    assert T.shape == C.shape == (toy.n_obs, len(exp.targets))
    # a cell is never both control and ring-exposed for the same target
    for R in rings:
        assert (R.multiply(C)).nnz == 0
    # controls under the clean policy have no perturbed neighbours
    r, _ = C.nonzero()
    assert (exp.any_perturbed_within[r] == 0).all()
    # autonomous-treated cells carry the target
    r, k = T.nonzero()
    assert all(exp.labels[i] == exp.targets[j] for i, j in zip(r, k))


def test_e1_runs_and_has_expected_rows(toy):
    exp = compute_exposure(toy, [0, 20, 40])
    Y = lognorm(toy)[:, :3]
    tab = E1Stratified(n_perm=5, min_cells=2, min_valid_perm=2).fit(toy, exp, Y, ["G0", "G1", "G2"])
    df = tab.df
    assert set(df["kind"]) == {"autonomous", "spillover"}
    assert df["ring"].max() == 1
    assert {"A", "B", "all"} <= set(df["cell_type"])
    ident = df[df["identified"]]
    assert len(ident) > 0
    assert ident["pvalue"].between(0, 1).all()
    q = tab.with_fdr()
    assert q.loc[q["identified"], "qvalue"].notna().all()


def test_e1_recovers_planted_autonomous_effect():
    a = make_toy(n=1500, n_genes=5, n_samples=2, seed=3, field_um=1500.0)
    # plant a strong autonomous effect of GENE_A on gene G0 (on the raw count scale)
    m = (a.obs["target"] == "GENE_A").to_numpy()
    a.X[m, 0] += 20
    exp = compute_exposure(a, [0, 15, 30])
    Y = lognorm(a)
    tab = E1Stratified(n_perm=20, min_cells=3).fit(a, exp, Y, list(a.var_names))
    df = tab.df
    row = df[
        (df.target == "GENE_A")
        & (df.kind == "autonomous")
        & (df.cell_type == "all")
        & (df.outcome == "G0")
    ].iloc[0]
    assert row["estimate"] > 0.5
    assert row["pvalue"] < 0.1
    assert row["pvalue_perm"] <= 1 / 21 + 1e-9  # floor with 20 permutations
    tab_p = E1Stratified(n_perm=20, min_cells=3, pvalue="perm").fit(a, exp, Y, list(a.var_names))
    assert (tab_p.df.loc[tab_p.df.identified, "pvalue"] >= 1 / 21 - 1e-9).all()
    # and no autonomous effect on G1 for GENE_A at the same scale
    row1 = df[
        (df.target == "GENE_A")
        & (df.kind == "autonomous")
        & (df.cell_type == "all")
        & (df.outcome == "G1")
    ].iloc[0]
    assert abs(row1["estimate"]) < row["estimate"]


def test_ntc_pseudo_targets_present(toy):
    exp = compute_exposure(toy, [0, 20])
    assert any(t.startswith("NTC:") for t in exp.targets)
    exp2 = compute_exposure(toy, [0, 20], include_ntc=False)
    assert not any(t.startswith("NTC:") for t in exp2.targets)


def test_exposure_recompute_roundtrip(toy):
    exp = compute_exposure(toy, [0, 20, 40])
    before = exp.total_count_matrix().toarray().copy()
    lab = exp.labels.copy()
    rng = np.random.default_rng(0)
    exp.recompute(rng.permutation(lab))
    exp.recompute(lab)
    assert np.array_equal(before, exp.total_count_matrix().toarray())
