"""Test of the undetected-sibling model on the simulator (pre-registered 2026-10-05 09:35).

Scenarios: clonal assignment and barcode misassignment, with a fraction of guide calls hidden
(detection efficiency 40%, 70%, 100%). For each spilling target the unassigned-neighbour
profile (uncorrected) and the sibling-corrected profile are compared with the planted
spillover profile (effective LFC per unit dose, sign and correlation) and with the planted
autonomous profile (contamination). Also reports whether pi tracks the hidden fraction.

Usage: uv run python scripts/benchmark_sibling.py [--n-cells 20000] [--reps 2]
Writes results/<hash>/sibling_benchmark.csv and .md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.config import results_dir
from spatialspill.exposure import compute_exposure
from spatialspill.loaders.registry import load_dataset
from spatialspill.outcomes import lognorm
from spatialspill.sibling import sibling_corrected_spillover
from spatialspill.simulator import (
    SimConfig,
    geometry_from,
    simulate,
    synthetic_geometry,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cells", type=int, default=20000)
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    cfg = {"script": "benchmark_sibling", "n_cells": a.n_cells, "reps": a.reps, "seed": a.seed}
    out = results_dir(cfg)
    try:
        geom = geometry_from(load_dataset("perturb_fish"), a.n_cells, seed=a.seed)
    except Exception:
        geom = synthetic_geometry(a.n_cells, seed=a.seed)
    scenarios = {
        "clonal": {"clonal_radius_um": 40.0, "clonal_size": 6},
        "misassign": {"p_misassign": 0.3},
        "clonal_misassign": {"clonal_radius_um": 40.0, "clonal_size": 6, "p_misassign": 0.3},
    }
    rows = []
    for sname, delta in scenarios.items():
        for det in (0.4, 0.7, 1.0):
            for rep in range(a.reps):
                sc = SimConfig(
                    n_genes=60,
                    n_targets=10,
                    frac_assigned=0.3,
                    frac_ntc_of_assigned=0.3,
                    frac_auto_nonzero=0.3,
                    auto_lfc_sd=1.5,
                    frac_spill_targets=1.0,
                    frac_spill_genes=0.2,
                    spill_lfc_sd=1.5,
                    spill_scale_um=30.0,
                    seed=a.seed * 100 + rep,
                    **delta,
                )
                sim = simulate(geom, sc)
                rng = np.random.default_rng(rep)
                assigned = (sim.obs["guide"] != "none").to_numpy()
                hide = assigned & (rng.random(sim.n_obs) >= det)
                sim.obs.loc[hide, ["guide", "target"]] = "none"
                sim.obs.loc[hide, "is_ntc"] = False
                sim.obs.loc[hide, "is_perturbed"] = False
                exp = compute_exposure(sim, [0, 15, 30, 60])
                Y = lognorm(sim)
                truth = sim.uns["truth"]
                p = np.asarray(truth["base_props"])
                for k, tgt in enumerate(truth["targets"]):
                    tau_s = np.asarray(truth["tau_spill"][k])
                    tau_a = np.asarray(truth["tau_auto"][k])
                    if not np.any(tau_s):
                        continue
                    eff_s = tau_s - np.log(np.sum(p * np.exp(tau_s)))
                    eff_a = tau_a - np.log(np.sum(p * np.exp(tau_a)))
                    est = sibling_corrected_spillover(
                        sim, exp, Y, tgt, ring=None, reference="ntc", n_boot=30, seed=rep
                    )
                    if est is None:
                        continue
                    nz = tau_s != 0
                    rows.append(
                        {
                            "scenario": sname,
                            "detection": det,
                            "rep": rep,
                            "target": tgt,
                            "pi": est.pi,
                            "n_unassigned_nb": est.n_unassigned_neighbours,
                            "r_unc_spill": float(np.corrcoef(est.uncorrected, eff_s)[0, 1]),
                            "r_cor_spill": float(np.corrcoef(est.corrected, eff_s)[0, 1]),
                            "r_unc_auto": float(np.corrcoef(est.uncorrected, eff_a)[0, 1]),
                            "r_cor_auto": float(np.corrcoef(est.corrected, eff_a)[0, 1]),
                            "sign_unc": float(
                                (np.sign(est.uncorrected[nz]) == np.sign(eff_s[nz])).mean()
                            ),
                            "sign_cor": float(
                                (np.sign(est.corrected[nz]) == np.sign(eff_s[nz])).mean()
                            ),
                        }
                    )
                print(sname, det, rep, len(rows), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "sibling_benchmark.csv", index=False)
    summ = (
        df.groupby(["scenario", "detection"])[
            ["pi", "r_unc_spill", "r_cor_spill", "r_unc_auto", "r_cor_auto", "sign_unc", "sign_cor"]
        ]
        .mean()
        .round(3)
    )
    (out / "sibling_benchmark.md").write_text(
        "# Sibling model benchmark\n\n" + summ.to_markdown() + "\n"
    )
    print(summ.to_string())
    print(f"results in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
