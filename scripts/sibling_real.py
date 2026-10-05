"""Exploratory: undetected-sibling model on a real dataset.

Usage: uv run python scripts/sibling_real.py --dataset perturb_fish --compare results/<e2 run>
Writes results/<hash>/sibling_real.csv and .json: pi_hat per target with bootstrap CI, number of
unassigned neighbours, and (if --compare is given) sign agreement and correlation between the
corrected profiles and the E2 ring-1/2 NTC-recipient estimates for the same targets.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.config import results_dir
from spatialspill.exposure import compute_exposure
from spatialspill.loaders.registry import load_dataset
from spatialspill.outcomes import lognorm
from spatialspill.sibling import sibling_table


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="perturb_fish")
    ap.add_argument("--compare", default=None)
    ap.add_argument("--bins", nargs="*", type=float, default=[0, 15, 30, 60])
    a = ap.parse_args()
    cfg = {"script": "sibling_real", "dataset": a.dataset, "bins": a.bins, "compare": a.compare}
    out = results_dir(cfg)
    adata = load_dataset(a.dataset)
    exp = compute_exposure(adata, a.bins)
    Y = lognorm(adata)
    tab = sibling_table(adata, exp, Y, [str(v) for v in adata.var_names], ring=None, n_boot=200)
    tab.to_csv(out / "sibling_real.csv", index=False)
    per = (
        tab.groupby("target")
        .agg(
            pi=("pi", "first"),
            pi_low=("pi_low", "first"),
            pi_high=("pi_high", "first"),
            n_unassigned_nb=("n_unassigned_neighbours", "first"),
        )
        .reset_index()
    )
    summ: dict[str, object] = {
        "n_targets": len(per),
        "pi_median": float(per["pi"].median()),
        "pi_iqr": [float(per["pi"].quantile(0.25)), float(per["pi"].quantile(0.75))],
        "pi_max": float(per["pi"].max()),
    }
    if a.compare:
        f = glob.glob(str(Path(a.compare) / "*_e2_estimates.csv"))
        if f:
            e2 = pd.read_csv(f[0])
            e2 = e2[
                (e2["kind"] == "spillover")
                & (e2["cell_type"] == "all")
                & e2["ring"].isin([1, 2])
                & e2["identified"]
            ]
            e2 = e2.groupby(["target", "outcome"], as_index=False)["estimate"].mean()
            m = tab.merge(e2, on=["target", "outcome"], how="inner")
            agree_cor = float((np.sign(m["corrected"]) == np.sign(m["estimate"])).mean())
            agree_unc = float((np.sign(m["uncorrected"]) == np.sign(m["estimate"])).mean())
            r_cor = float(np.corrcoef(m["corrected"], m["estimate"])[0, 1])
            r_unc = float(np.corrcoef(m["uncorrected"], m["estimate"])[0, 1])
            summ.update(
                {
                    "n_pairs_compared": len(m),
                    "sign_agreement_corrected": agree_cor,
                    "sign_agreement_uncorrected": agree_unc,
                    "r_corrected": r_cor,
                    "r_uncorrected": r_unc,
                }
            )
    (out / "sibling_real.json").write_text(json.dumps(summ, indent=2) + "\n")
    per.to_csv(out / "sibling_pi_per_target.csv", index=False)
    print(json.dumps(summ, indent=2))
    print(per.sort_values("pi", ascending=False).head(12).to_string())
    print(f"results in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
