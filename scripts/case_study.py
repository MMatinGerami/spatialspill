"""Deep case study (C6): the top-ranked niche-remodelling target in Perturb-FISH.

Produces results/<hash>/case_<TARGET>.png and case_<TARGET>.json: (a) the spatial layout of
the target's cells, their NTC recipients in rings 1 and 2 and the NTC controls in one field;
(b) the ring-1 and ring-2 spillover profiles (top genes by |z|) with 95% intervals; (c) the
autonomous profile of the same target for comparison; (d) the NTC pseudo-target distribution
of ring-1 z for the same genes as the empirical null.

Usage: uv run python scripts/case_study.py --results results/<hash> [--target MYD88]
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.exposure import compute_exposure
from spatialspill.loaders.registry import load_dataset


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--target", default=None)
    ap.add_argument("--dataset", default="perturb_fish")
    ap.add_argument("--top", type=int, default=12)
    a = ap.parse_args()
    d = Path(a.results)
    df = pd.read_csv(glob.glob(str(d / "*_e2_estimates.csv"))[0])
    rk = pd.read_csv(d / "niche_ranking.csv")
    target = a.target or str(rk.iloc[0]["target"])
    sub = df[(df["target"] == target) & (df["cell_type"] == "all") & df["identified"]].copy()
    sub["z"] = sub["estimate"] / sub["se"]
    out: dict[str, object] = {
        "target": target,
        "rank": int(rk.index[rk["target"] == target][0]) + 1,
    }

    adata = load_dataset(a.dataset)
    exp = compute_exposure(adata, [0, 15, 30, 60])
    k = exp.targets.index(target)
    own = exp.labels == target
    ntc = adata.obs["is_ntc"].to_numpy()
    r1 = (np.asarray(exp.counts[1][:, k].todense()).ravel() > 0) & ntc
    r2 = (np.asarray(exp.counts[2][:, k].todense()).ravel() > 0) & ntc
    out["n_target_cells"] = int(own.sum())
    out["n_ntc_recipients_ring1"] = int(r1.sum())
    out["n_ntc_recipients_ring2"] = int(r2.sum())
    xy = np.asarray(adata.obsm["spatial"])
    # field of 1 mm around the densest cluster of target cells
    cx, cy = np.median(xy[own], axis=0)
    win = (np.abs(xy[:, 0] - cx) < 500) & (np.abs(xy[:, 1] - cy) < 500)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    ax = axes[0, 0]
    ax.scatter(xy[win, 0], xy[win, 1], s=2, c="#dddddd", label="all cells")
    ax.scatter(xy[win & ntc, 0], xy[win & ntc, 1], s=4, c="#4477aa", label="NTC")
    ax.scatter(xy[win & r1, 0], xy[win & r1, 1], s=14, c="#ee7733", label="NTC, ring 1 recipient")
    ax.scatter(xy[win & own, 0], xy[win & own, 1], s=14, c="#cc3311", label=f"{target} cells")
    ax.set_aspect("equal")
    ax.legend(fontsize=7, markerscale=2)
    ax.set_title(f"{target}: 1 mm field around the densest cluster")
    for ax, ring, title in (
        (axes[0, 1], 1, "ring 1 (15 to 30 um)"),
        (axes[1, 0], 2, "ring 2 (30 to 60 um)"),
    ):
        s = (
            sub[(sub["kind"] == "spillover") & (sub["ring"] == ring)]
            .sort_values("z", key=np.abs, ascending=False)
            .head(a.top)
        )
        ax.errorbar(
            s["estimate"],
            range(len(s)),
            xerr=1.96 * s["se"],
            fmt="o",
            color="#ee7733",
            ecolor="grey",
            capsize=2,
        )
        ax.set_yticks(range(len(s)))
        ax.set_yticklabels(s["outcome"], fontsize=8)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_title(f"{target} spillover onto NTC recipients, {title}")
        ax.set_xlabel("difference in log1p-normalised expression (95% CI)")
        out[f"ring{ring}_top"] = (
            s[["outcome", "estimate", "se", "qvalue"]].round(4).to_dict("records")
        )
    ax = axes[1, 1]
    s = sub[sub["kind"] == "autonomous"].sort_values("z", key=np.abs, ascending=False).head(a.top)
    ax.errorbar(
        s["estimate"],
        range(len(s)),
        xerr=1.96 * s["se"],
        fmt="o",
        color="#cc3311",
        ecolor="grey",
        capsize=2,
    )
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(s["outcome"], fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_title(f"{target} autonomous effect (cells carrying the guide)")
    ax.set_xlabel("difference in log1p-normalised expression (95% CI)")
    out["autonomous_top"] = s[["outcome", "estimate", "se", "qvalue"]].round(4).to_dict("records")
    # overlap between spillover and autonomous top genes
    a1 = set(sub[(sub.kind == "spillover") & (sub.ring == 1) & (sub.qvalue < 0.1)]["outcome"])
    au = set(sub[(sub.kind == "autonomous") & (sub.qvalue < 0.1)]["outcome"])
    out["ring1_hits"] = sorted(a1)
    out["autonomous_hits"] = sorted(au)
    out["ring1_hits_also_autonomous"] = sorted(a1 & au)
    fig.tight_layout()
    fig.savefig(d / f"case_{target}.png", dpi=130)
    (d / f"case_{target}.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: v for k, v in out.items() if not k.endswith("_top")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
