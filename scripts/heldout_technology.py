"""Held-out-technology test (C5): train a predictor of a target's niche-remodelling score on one
dataset and evaluate it on another technology, through a shared gene embedding (GO, symbols
upper-cased so that human and mouse orthologues align).

Usage: uv run python scripts/heldout_technology.py --train results/<fish> --test results/<multi>
Writes results/summary/heldout_technology.json. The score is Spearman correlation between
predicted and observed niche scores on the test dataset, against a permuted-embedding null.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.kernel_ridge import KernelRidge

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.embeddings import go_embedding

SUPPORT = Path(__file__).resolve().parents[1] / "data" / "raw" / "supporting"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--n-null", type=int, default=200)
    a = ap.parse_args()
    tr = pd.read_csv(Path(a.train) / "niche_ranking.csv")
    te = pd.read_csv(Path(a.test) / "niche_ranking.csv")
    E = pd.concat(
        [
            go_embedding(SUPPORT / "go" / "goa_human.gaf.gz", SUPPORT / "go" / "go-basic.obo"),
            go_embedding(SUPPORT / "go" / "goa_mouse.gaf.gz", SUPPORT / "go" / "go-basic.obo"),
        ]
    )
    E = E[~E.index.duplicated()]
    tr["t"] = tr["target"].str.upper()
    te["t"] = te["target"].str.upper()
    tr = tr[tr["t"].isin(E.index)]
    te = te[te["t"].isin(E.index)]
    overlap = sorted(set(tr["t"]) & set(te["t"]))
    Xtr = E.loc[tr["t"]].to_numpy()
    Xte = E.loc[te["t"]].to_numpy()
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
    ytr = tr["score"].to_numpy()
    yte = te["score"].to_numpy()
    model = KernelRidge(alpha=1.0, kernel="rbf", gamma=1.0 / Xtr.shape[1]).fit(Xtr, ytr)
    pred = model.predict(Xte)
    rho = stats.spearmanr(pred, yte).statistic
    rng = np.random.default_rng(0)
    null = []
    for _ in range(a.n_null):
        perm = rng.permutation(len(Xtr))
        m = KernelRidge(alpha=1.0, kernel="rbf", gamma=1.0 / Xtr.shape[1]).fit(Xtr[perm], ytr)
        null.append(stats.spearmanr(m.predict(Xte), yte).statistic)
    null_arr = np.array(null, dtype=float)
    out = {
        "train": str(a.train),
        "test": str(a.test),
        "n_train_targets": len(tr),
        "n_test_targets": len(te),
        "n_overlapping_targets": len(overlap),
        "overlapping_targets": overlap,
        "spearman": float(rho),
        "null_mean": float(np.nanmean(null_arr)),
        "null_sd": float(np.nanstd(null_arr)),
        "empirical_p": float((np.sum(null_arr >= rho) + 1) / (len(null_arr) + 1)),
    }
    Path("results/summary").mkdir(parents=True, exist_ok=True)
    Path("results/summary/heldout_technology.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
