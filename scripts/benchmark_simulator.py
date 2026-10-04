"""Benchmark estimators on simulator scenarios (C2). Writes results/<hash>/benchmark_*.csv.

Usage: uv run python scripts/benchmark_simulator.py --config configs/benchmark_small.yaml
Geometry comes from a real dataset when available (Perturb-FISH tumour subsample), else synthetic.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.benchmark import score
from spatialspill.config import load_config, results_dir
from spatialspill.estimators import E2GLM, E1Stratified, E3DoublyRobust
from spatialspill.estimators.groups import GroupConfig
from spatialspill.exposure import compute_exposure
from spatialspill.outcomes import lognorm
from spatialspill.simulator import (
    SimConfig,
    default_scenarios,
    geometry_from,
    simulate,
    synthetic_geometry,
)


def get_geometry(cfg) -> tuple:
    src = cfg.get("geometry", "synthetic")
    if src == "synthetic":
        return synthetic_geometry(int(cfg.n_cells), seed=int(cfg.seed)), "synthetic"
    from spatialspill.loaders.registry import load_dataset

    try:
        a = load_dataset(src)
        return geometry_from(a, int(cfg.n_cells), seed=int(cfg.seed)), src
    except Exception as exc:
        print(f"geometry {src} unavailable ({exc}); using synthetic", file=sys.stderr)
        return synthetic_geometry(int(cfg.n_cells), seed=int(cfg.seed)), "synthetic"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("overrides", nargs="*")
    a = ap.parse_args()
    cfg = load_config(a.config, a.overrides)
    out = results_dir(cfg)
    geom, geom_name = get_geometry(cfg)
    scen = {s.name: s for s in default_scenarios(int(cfg.seed))}
    names = list(cfg.get("scenarios", list(scen)))
    all_scores, meta = [], []
    for name in names:
        s = scen[name]
        overrides = {
            "n_genes": int(cfg.n_genes),
            "n_targets": int(cfg.n_targets),
            "frac_assigned": float(cfg.frac_assigned),
        }
        overrides.update(dict(cfg.get("sim_overrides", {})))
        sc = SimConfig(**{**s.cfg.__dict__, **overrides})
        for rep in range(int(cfg.n_reps)):
            sc.seed = int(cfg.seed) * 1000 + rep
            t0 = time.time()
            sim = simulate(geom, sc)
            exp = compute_exposure(sim, list(cfg.distance_bins_um))
            Y = lognorm(sim)
            for est_name in cfg.estimators:
                gcfg = GroupConfig(
                    control_policy=str(cfg.control_policy),
                    clean_controls=bool(cfg.clean_controls),
                )
                if est_name == "E1":
                    est = E1Stratified(
                        n_perm=int(cfg.n_perm),
                        seed=sc.seed,
                        min_cells=int(cfg.min_cells),
                        groups=gcfg,
                    )
                elif est_name == "E1_analytic":
                    est = E1Stratified(
                        n_perm=int(cfg.n_perm),
                        seed=sc.seed,
                        min_cells=int(cfg.min_cells),
                        groups=gcfg,
                        ci="analytic",
                    )
                elif est_name == "E2":
                    est = E2GLM(min_cells=int(cfg.min_cells), groups=gcfg)
                elif est_name == "E3":
                    est = E3DoublyRobust(
                        n_draws=int(cfg.n_perm),
                        seed=sc.seed,
                        min_cells=int(cfg.min_cells),
                        groups=gcfg,
                        n_boot=100,
                    )
                else:
                    raise KeyError(est_name)
                tab = est.fit(sim, exp, Y, list(sim.var_names)).with_fdr()
                sc_df = score(tab, sim, float(cfg.fdr))
                sc_df.insert(0, "estimator", est_name)
                sc_df.insert(0, "rep", rep)
                sc_df.insert(0, "scenario", name)
                all_scores.append(sc_df)
            meta.append(
                {
                    "scenario": name,
                    "rep": rep,
                    "seconds": round(time.time() - t0, 1),
                    "n_cells": sim.n_obs,
                    "n_perturbed": int(sim.obs.is_perturbed.sum()),
                    "n_ntc": int(sim.obs.is_ntc.sum()),
                    "n_misassigned": int(sim.uns["truth"]["n_misassigned"]),
                }
            )
            print(json.dumps(meta[-1]), flush=True)
    res = pd.concat(all_scores, ignore_index=True)
    res.to_csv(out / "benchmark_scores.csv", index=False)
    pd.DataFrame(meta).to_csv(out / "benchmark_meta.csv", index=False)
    summary = (
        res.groupby(["scenario", "estimator", "kind", "ring"])
        .agg(
            null_fpr_q=("null_fpr_q", "mean"),
            null_cov95=("null_coverage95", "mean"),
            ntc_cov95=("ntc_coverage95", "mean"),
            ntc_fpr_p05=("ntc_fpr_p05", "mean"),
            power_q=("power_q", "mean"),
            auroc=("auroc", "mean"),
            fdp_q=("fdp_q", "mean"),
            n_tests=("n_tests", "mean"),
        )
        .round(3)
        .reset_index()
    )
    summary.to_csv(out / "benchmark_summary.csv", index=False)
    (out / "benchmark_summary.md").write_text(
        f"# Simulator benchmark (geometry: {geom_name})\n\n"
        + summary.to_markdown(index=False)
        + "\n"
    )
    print(summary.to_string())
    print(f"results in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
