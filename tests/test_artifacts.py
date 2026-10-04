import numpy as np
import pandas as pd

from spatialspill.artifacts import artifact_fraction, bleedthrough_table, covariate_shift


def _table(alpha0: float, alpha2: float, n_out: int = 30, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for t in ["T1", "T2", "T3"]:
        a = rng.normal(0, 1, n_out)
        for g in range(n_out):
            rows.append(
                {
                    "target": t,
                    "kind": "autonomous",
                    "ring": -1,
                    "cell_type": "all",
                    "outcome": f"g{g}",
                    "estimate": a[g],
                    "identified": True,
                    "qvalue": 0.5,
                }
            )
        for ring, alpha in [(0, alpha0), (1, 0.5 * alpha0), (2, alpha2)]:
            s = alpha * a + rng.normal(0, 0.1, n_out)
            for g in range(n_out):
                rows.append(
                    {
                        "target": t,
                        "kind": "spillover",
                        "ring": ring,
                        "cell_type": "all",
                        "outcome": f"g{g}",
                        "estimate": s[g],
                        "identified": True,
                        "qvalue": 0.5,
                    }
                )
    return pd.DataFrame(rows)


def test_bleedthrough_recovers_alpha_and_fraction():
    df = _table(alpha0=0.3, alpha2=0.0)
    bt = bleedthrough_table(df)
    r0 = bt[bt.ring == 0]
    assert np.allclose(r0["alpha_hat"], 0.3, atol=0.05)
    assert (r0["r2"] > 0.8).all()
    summ = artifact_fraction(bt)
    assert summ["r2_ring0_minus_last"] > 0.6
    assert summ["alpha_ring0"] > 0.25 and abs(summ["alpha_last"]) < 0.1


def test_no_bleedthrough_gives_small_fraction():
    df = _table(alpha0=0.0, alpha2=0.0)
    summ = artifact_fraction(bleedthrough_table(df))
    assert abs(summ["r2_ring0_minus_last"]) < 0.3


def test_empty_and_covariate_shift():
    assert artifact_fraction(pd.DataFrame())["n_targets"] == 0
    df = _table(0.2, 0.0)
    shift = covariate_shift(df, df)
    assert (shift["mean_abs_change"] == 0).all() and (shift["corr"] > 0.99).all()
