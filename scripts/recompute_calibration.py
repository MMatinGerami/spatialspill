"""Recompute <dataset>_calibration.json for every results directory from its estimate table
(used after the summary definition changed; the estimates themselves are untouched)."""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_pipeline import calibration_summary

for d in sorted(Path("results").glob("*/")):
    cfg = d / "config.yaml"
    if not cfg.exists():
        continue
    fdr = float(yaml.safe_load(cfg.read_text()).get("fdr", 0.1))
    for f in glob.glob(str(d / "*_e*_estimates.csv")):
        name = Path(f).name.split("_e")[0]
        df = pd.read_csv(f)
        summ = calibration_summary(df, fdr)
        summ["dataset"] = name
        summ["estimator"] = Path(f).name.split("_")[-2].upper()
        (d / f"{name}_calibration.json").write_text(json.dumps(summ, indent=2) + "\n")
        print(
            d.name, name, summ["estimator"], summ["n_tests"], summ.get("ntc_fraction_p_below_0.05")
        )
