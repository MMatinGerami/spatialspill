"""Power analysis (pre-registered 2026-10-05 09:35): power to detect planted spillover as a
function of the number of confirmed non-target recipients per target and of effect size.

Usage: uv run python scripts/power_analysis.py --config configs/power_grid.yaml
Writes results/<hash>/power_grid.csv and power_grid.png (power at q<0.1 in rings 1 and 2 vs
median NTC recipients per target, one panel per effect size; NTC false positives as a check).
"""

from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.benchmark import score
from spatialspill.config import load_config, results_dir
from spatialspill.estimators import E2GLM, E1Stratified
from spatialspill.estimators.groups import GroupConfig
from spatialspill.exposure import compute_exposure
from spatialspill.loaders.registry import load_dataset
from spatialspill.outcomes import lognorm
from spatialspill.simulator import (
    SimConfig,
    geometry_from,
    simulate,
    synthetic_geometry,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("overrides", nargs="*")
    a = ap.parse_args()
    cfg = load_config(a.config, a.overrides)
    out = results_dir(cfg)
    try:
        geom = geometry_from(load_dataset(str(cfg.geometry)), int(cfg.n_cells), seed=int(cfg.seed))
    except Exception as exc:
        print(f"geometry unavailable ({exc}); synthetic", file=sys.stderr)
        geom = synthetic_geometry(int(cfg.n_cells), seed=int(cfg.seed))
    gcfg = GroupConfig(control_policy="ntc", clean_controls=bool(cfg.clean_controls))
    rows = []
    grid = list(
        itertools.product(list(cfg.frac_assigned), list(cfg.frac_ntc), list(cfg.spill_lfc_sd))
    )
    for i, (fa, fn, sd) in enumerate(grid):
        sc = SimConfig(
            n_genes=int(cfg.n_genes),
            n_targets=int(cfg.n_targets),
            frac_assigned=float(fa),
            frac_ntc_of_assigned=float(fn),
            spill_lfc_sd=float(sd),
            frac_spill_genes=0.2,
            frac_spill_targets=0.5,
            spill_scale_um=30.0,
            seed=int(cfg.seed) * 100 + i,
        )
        t0 = time.time()
        sim = simulate(geom, sc)
        tile = float(cfg.tile_um)
        xy = np.asarray(sim.obsm["spatial"], dtype=float)
        sim.obs["tile"] = [
            f"{int(x)}_{int(y)}"
            for x, y in zip(np.floor(xy[:, 0] / tile), np.floor(xy[:, 1] / tile))
        ]
        exp = compute_exposure(sim, list(cfg.distance_bins_um))
        Y = lognorm(sim)
        ntc = sim.obs["is_ntc"].to_numpy()
        tot = exp.total_count_matrix().toarray()
        rec = [
            (tot[ntc, k] > 0).sum() for k, t in enumerate(exp.targets) if not t.startswith("NTC:")
        ]
        med_rec = float(np.median(rec))
        ests = {
            "E1_tile": E1Stratified(
                n_perm=int(cfg.n_perm),
                seed=sc.seed,
                min_cells=5,
                groups=gcfg,
                strata_keys=("sample", "cell_type", "tile"),
            ),
            "E2_spatial_perm": E2GLM(
                min_cells=5,
                groups=gcfg,
                spatial_basis=int(cfg.spatial_basis),
                n_perm=int(cfg.n_perm),
                seed=sc.seed,
            ),
        }
        for name, est in ests.items():
            tab = est.fit(sim, exp, Y, list(sim.var_names)).with_fdr()
            sc_df = score(tab, sim, 0.1)
            for _, r in sc_df.iterrows():
                rows.append(
                    {
                        "frac_assigned": fa,
                        "frac_ntc": fn,
                        "spill_lfc_sd": sd,
                        "median_ntc_recipients": med_rec,
                        "n_ntc": int(ntc.sum()),
                        "estimator": name,
                        **r.to_dict(),
                    }
                )
        print(
            f"{i + 1}/{len(grid)} fa={fa} fn={fn} sd={sd} ntc_rec={med_rec:.0f} {time.time() - t0:.0f}s",
            flush=True,
        )
        pd.DataFrame(rows).to_csv(out / "power_grid.csv", index=False)
    df = pd.DataFrame(rows)
    sp = df[(df["kind"] == "spillover") & (df["ring"] > 0)]
    sds = sorted(sp["spill_lfc_sd"].unique())
    fig, axes = plt.subplots(1, len(sds), figsize=(5 * len(sds), 4.2), sharey=True)
    for ax, sd in zip(np.atleast_1d(axes), sds):
        for name, g in sp[sp["spill_lfc_sd"] == sd].groupby("estimator"):
            gg = (
                g.groupby("median_ntc_recipients", as_index=False)["power_q"]
                .mean()
                .sort_values("median_ntc_recipients")
            )
            ax.plot(gg["median_ntc_recipients"], gg["power_q"], "o-", label=name)
        ax.set_title(f"planted spillover LFC sd {sd}")
        ax.set_xlabel("median confirmed NTC recipients per target")
        ax.set_xscale("log")
        ax.axhline(0.8, color="k", ls="--", lw=0.8)
    np.atleast_1d(axes)[0].set_ylabel("power at q < 0.1 (rings 1 and 2)")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "power_grid.png", dpi=140)
    summ = (
        sp.groupby(["estimator", "spill_lfc_sd", "frac_assigned", "frac_ntc"])
        .agg(
            ntc_recipients=("median_ntc_recipients", "first"),
            power=("power_q", "mean"),
            ntc_fpr=("ntc_fpr_p05", "mean"),
            fdp=("fdp_q", "mean"),
        )
        .round(3)
        .reset_index()
    )
    summ.to_csv(out / "power_summary.csv", index=False)
    print(summ.to_string())
    print(f"results in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
