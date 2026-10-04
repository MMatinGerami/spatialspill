"""Run the estimator pipeline for the datasets in a config and write results/<hash>/.

Usage: uv run python scripts/run_pipeline.py --config configs/smoke.yaml [key=value ...]

Outputs per dataset:
  <dataset>_e1_estimates.csv   full E1 table with BH q-values
  <dataset>_calibration.json   false-positive rate of NTC pseudo-targets at the configured FDR,
                               plus counts of identified estimands
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.config import load_config, results_dir
from spatialspill.estimators import E1Stratified
from spatialspill.estimators.groups import GroupConfig
from spatialspill.exposure import compute_exposure
from spatialspill.loaders.registry import load_dataset
from spatialspill.outcomes import lognorm


def calibration_summary(df, fdr: float) -> dict[str, object]:
    ident = df[df["identified"]]
    ntc = ident[ident["target"].str.startswith("NTC:")]
    tgt = ident[~ident["target"].str.startswith("NTC:")]
    out: dict[str, object] = {
        "n_tests": len(ident),
        "n_ntc_tests": len(ntc),
        "n_target_tests": len(tgt),
        "fdr": fdr,
        "ntc_fraction_q_below_fdr": float((ntc["qvalue"] < fdr).mean()) if len(ntc) else None,
        "ntc_fraction_p_below_0.05": float((ntc["pvalue"] < 0.05).mean()) if len(ntc) else None,
        "target_fraction_q_below_fdr": float((tgt["qvalue"] < fdr).mean()) if len(tgt) else None,
    }
    for kind in ("autonomous", "spillover"):
        sub = tgt[tgt["kind"] == kind]
        out[f"n_{kind}_hits"] = int((sub["qvalue"] < fdr).sum())
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("overrides", nargs="*")
    a = ap.parse_args()
    cfg = load_config(a.config, a.overrides)
    out = results_dir(cfg)
    rng_seed = int(cfg.seed)
    for name in cfg.datasets:
        adata = load_dataset(name)
        print(f"[{name}] {adata.n_obs} units x {adata.n_vars} features", flush=True)
        bins = list(cfg.distance_bins_um)
        exp = compute_exposure(adata, bins)
        Y = lognorm(adata)
        if "max_outcomes" in cfg and cfg.max_outcomes and Y.shape[1] > int(cfg.max_outcomes):
            tot = np.asarray(adata.X.sum(axis=0)).ravel()
            keep = np.sort(np.argsort(-tot)[: int(cfg.max_outcomes)])
            Y = Y[:, keep]
            names = [str(v) for v in adata.var_names[keep]]
        else:
            names = [str(v) for v in adata.var_names]
        est = E1Stratified(
            n_perm=int(cfg.n_perm),
            seed=rng_seed,
            strata_keys=tuple(cfg.strata),
            min_cells=int(cfg.get("min_cells", 5)),
            groups=GroupConfig(
                control_policy=str(cfg.get("control_policy", "ntc")),
                clean_controls=bool(cfg.get("clean_controls", True)),
            ),
        )
        tab = est.fit(adata, exp, Y, names)
        df = tab.with_fdr()
        df.to_csv(out / f"{name}_e1_estimates.csv", index=False)
        summ = calibration_summary(df, float(cfg.fdr))
        summ["dataset"] = name
        summ["targets"] = len(exp.targets)
        (out / f"{name}_calibration.json").write_text(json.dumps(summ, indent=2) + "\n")
        print(json.dumps(summ, indent=2))
    print(f"results in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
