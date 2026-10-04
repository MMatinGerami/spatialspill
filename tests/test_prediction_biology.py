import numpy as np
import pandas as pd

from spatialspill.biology import load_lr_pairs, lr_enrichment, rank_niche_genes
from spatialspill.prediction import loto_predict, response_matrix


def _est(seed=0, n_t=12, n_o=15):
    rng = np.random.default_rng(seed)
    rows = []
    for t in range(n_t):
        for o in range(n_o):
            for kind, ring in (("autonomous", -1), ("spillover", 1)):
                e = rng.normal(0, 1)
                rows.append(
                    {
                        "target": f"T{t}",
                        "kind": kind,
                        "ring": ring,
                        "cell_type": "all",
                        "outcome": f"O{o}",
                        "estimate": e,
                        "se": 0.5,
                        "pvalue": 0.5,
                        "qvalue": 0.5 if abs(e) < 2 else 0.01,
                        "identified": True,
                    }
                )
    rows.append(
        {
            "target": "NTC:g1",
            "kind": "spillover",
            "ring": 1,
            "cell_type": "all",
            "outcome": "O0",
            "estimate": 0.1,
            "se": 0.5,
            "pvalue": 0.8,
            "qvalue": 0.9,
            "identified": True,
        }
    )
    return pd.DataFrame(rows)


def test_response_and_loto():
    df = _est()
    R = response_matrix(df, "spillover", rings=[1])
    assert R.shape == (12, 15)
    # embedding that encodes the response exactly should predict well; random should not
    E_good = R.copy()
    E_good.columns = [f"e{i}" for i in range(E_good.shape[1])]
    good = loto_predict(R, E_good, n_null=3)
    rng = np.random.default_rng(1)
    E_bad = pd.DataFrame(rng.normal(size=(12, 8)), index=R.index)
    bad = loto_predict(R, E_bad, n_null=3)
    assert good["mean_r"] > bad["mean_r"]
    assert good["n_targets"] == 12


def test_lr_and_ranking(tmp_path):
    p = tmp_path / "cpdb.csv"
    pd.DataFrame({"partner_a": ["T0", "T1"], "partner_b": ["O0", "O1"]}).to_csv(p, index=False)
    lr = load_lr_pairs(cellphonedb_csv=p)
    assert ("T0", "O0") in lr
    df = _est()
    res = lr_enrichment(df[df.kind == "spillover"], lr)
    assert (res["n_lr_pairs"] == 2 and 0 <= res["pvalue"] <= 1) or np.isnan(res["pvalue"])
    rk = rank_niche_genes(df, df)
    assert len(rk) == 12 and "score" in rk and "ntc_score_q95" in rk
