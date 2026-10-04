"""Per-dataset audit report (Phase 1): guide-call rate, cells per target, NTC counts, local
density, segmentation area, batch structure, spatial autocorrelation of guide labels (join
count against the stratified permutation null) and graph summaries for the default and
alternative neighbour graphs.

Usage: uv run python scripts/audit_datasets.py --dataset bundled [--out reports/audit]
Writes reports/audit/<dataset>/{audit.json, audit.md, *.png}.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.exposure import local_density
from spatialspill.graphs import GraphConfig, build_graph, degree
from spatialspill.loaders.registry import load_dataset
from spatialspill.permutation import permute_within_strata, strata_codes

DATASET_KWARGS: dict[str, dict] = {
    "perturb_multi": {"batches": ["4", "5", "8", "9", "10"], "max_cells_per_batch": 150_000},
}

# neighbourhood scale per dataset: (density radius, Delaunay max edge, radius-graph radius, knn)
# single-cell imaging: 30 um radius; Visium: spots sit on a 100 um hex lattice; DBiT: pixel units
SCALES: dict[str, tuple[float, float, float, int]] = {
    "default_um": (30.0, 50.0, 30.0, 10),
    "perturb_map": (110.0, 120.0, 110.0, 6),
    "pixel": (1.5, 1.5, 1.5, 8),
}


def join_count(A, labels: np.ndarray) -> float:
    """Number of graph edges joining two cells with the same (perturbed) label."""
    r, c = A.nonzero()
    m = r < c
    r, c = r[m], c[m]
    same = (labels[r] == labels[c]) & (labels[r] != "none") & (labels[r] != "NTC")
    return float(same.sum())


def audit(name: str, out_root: Path, n_perm: int = 100, seed: int = 0) -> dict:
    adata = load_dataset(name, **DATASET_KWARGS.get(name, {}))
    obs = adata.obs
    out = out_root / name
    out.mkdir(parents=True, exist_ok=True)
    units = adata.uns["spatialspill"].get("units", "um")
    rep: dict[str, object] = {
        "dataset": name,
        "technology": adata.uns["spatialspill"].get("technology"),
        "units": units,
        "n_units": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "n_samples": int(obs["sample"].nunique()),
        "guide_call_rate": float((obs["guide"] != "none").mean()),
        "n_perturbed": int(obs["is_perturbed"].sum()),
        "n_ntc": int(obs["is_ntc"].sum()),
        "n_ntc_guides": int(obs.loc[obs["is_ntc"], "guide"].nunique()),
        "n_targets": int(obs.loc[obs["is_perturbed"], "target"].nunique()),
        "cells_per_target": obs.loc[obs["is_perturbed"], "target"]
        .value_counts()
        .describe()
        .round(1)
        .to_dict(),
        "cell_types": obs["cell_type"].value_counts().to_dict(),
        "area_um2": obs["area"].describe().round(1).to_dict()
        if obs["area"].notna().any()
        else None,
        "total_counts": obs["total_counts"].describe().round(1).to_dict()
        if "total_counts" in obs
        else None,
    }
    if "n_guides" in obs:
        rep["guides_per_unit"] = obs["n_guides"].value_counts().sort_index().head(8).to_dict()
    if "n_barcodes_lenient" in obs:
        rep["barcodes_per_cell_lenient"] = (
            obs["n_barcodes_lenient"].value_counts().sort_index().head(8).to_dict()
        )

    # density
    scale_key = name if name in SCALES else ("pixel" if units != "um" else "default_um")
    radius, delaunay_max, radius_graph, knn = SCALES[scale_key]
    dens = local_density(adata, radius)
    rep["local_density"] = {"radius": radius, **pd.Series(dens).describe().round(2).to_dict()}
    obs["local_density"] = dens
    # density balance: perturbed vs NTC vs none, within strata
    rep["density_by_group"] = {
        "perturbed": float(dens[obs["is_perturbed"].to_numpy()].mean())
        if rep["n_perturbed"]
        else None,
        "ntc": float(dens[obs["is_ntc"].to_numpy()].mean()) if rep["n_ntc"] else None,
        "unassigned": float(dens[(obs["guide"] == "none").to_numpy()].mean()),
    }

    # graphs
    graphs = {}
    configs = {
        "delaunay_pruned": GraphConfig(kind="delaunay", max_edge_um=delaunay_max),
        "radius": GraphConfig(kind="radius", radius_um=radius_graph),
        f"knn{knn}": GraphConfig(kind="knn", n_neighs=knn),
    }
    strata = strata_codes(obs, ("sample", "cell_type"))
    labels = obs["target"].astype(str).to_numpy()
    rng = np.random.default_rng(seed)
    for gname, gcfg in configs.items():
        build_graph(adata, gcfg)
        A = adata.obsp["spatial_connectivities"]
        D = adata.obsp["spatial_distances"]
        deg = degree(adata)
        jc_obs = join_count(A, labels)
        null = np.array(
            [join_count(A, permute_within_strata(labels, strata, rng)) for _ in range(n_perm)]
        )
        graphs[gname] = {
            "mean_degree": float(deg.mean()),
            "isolated_fraction": float((deg == 0).mean()),
            "median_edge_length": float(np.median(D.data)) if D.nnz else None,
            "same_target_edges_observed": jc_obs,
            "same_target_edges_null_mean": float(null.mean()),
            "same_target_edges_null_sd": float(null.std()),
            "clonality_z": float((jc_obs - null.mean()) / null.std()) if null.std() > 0 else None,
        }
    rep["graphs"] = graphs
    rep["batch_structure"] = pd.crosstab(obs["sample"], obs["cell_type"]).to_dict()

    # figures
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    s0 = obs["sample"].iloc[0]
    m = (obs["sample"] == s0).to_numpy()
    xy = adata.obsm["spatial"][m]
    axes[0].scatter(xy[:, 0], xy[:, 1], s=1, c="lightgrey")
    pm = m & obs["is_perturbed"].to_numpy()
    axes[0].scatter(
        adata.obsm["spatial"][pm, 0],
        adata.obsm["spatial"][pm, 1],
        s=2,
        c="crimson",
        label="perturbed",
    )
    nm = m & obs["is_ntc"].to_numpy()
    axes[0].scatter(
        adata.obsm["spatial"][nm, 0], adata.obsm["spatial"][nm, 1], s=2, c="navy", label="NTC"
    )
    axes[0].set_title(f"{name}: sample {s0}")
    axes[0].set_aspect("equal")
    axes[0].legend(markerscale=4, fontsize=8)
    vc = obs.loc[obs["is_perturbed"], "target"].value_counts()
    axes[1].hist(vc.to_numpy(), bins=30, color="grey")
    axes[1].set_xlabel("units per target")
    axes[1].set_title(f"{len(vc)} targets")
    axes[2].hist(dens, bins=40, color="grey")
    axes[2].set_xlabel(f"neighbours within {radius} {units}")
    axes[2].set_title("local density")
    fig.tight_layout()
    fig.savefig(out / "overview.png", dpi=120)
    plt.close(fig)

    (out / "audit.json").write_text(json.dumps(rep, indent=2, default=str) + "\n")
    md = [
        f"# Audit: {name}",
        "",
        f"Generated by scripts/audit_datasets.py ({rep['technology']}, units: {units}).",
        "",
    ]
    md.append("| metric | value |\n|---|---|")
    for k in (
        "n_units",
        "n_features",
        "n_samples",
        "guide_call_rate",
        "n_perturbed",
        "n_ntc",
        "n_ntc_guides",
        "n_targets",
    ):
        v = rep[k]
        md.append(f"| {k} | {v:.4f} |" if isinstance(v, float) else f"| {k} | {v} |")
    md.append("")
    md.append("## Neighbour graphs and clonality\n")
    md.append(
        "| graph | mean degree | isolated | median edge | same-target edges obs | null mean (sd) | z |\n|---|---|---|---|---|---|---|"
    )
    for g, v in graphs.items():
        z = f"{v['clonality_z']:.1f}" if v["clonality_z"] is not None else "n/a"
        med = f"{v['median_edge_length']:.1f}" if v["median_edge_length"] is not None else "n/a"
        md.append(
            f"| {g} | {v['mean_degree']:.2f} | {v['isolated_fraction']:.3f} | {med} | {v['same_target_edges_observed']:.0f} | {v['same_target_edges_null_mean']:.1f} ({v['same_target_edges_null_sd']:.1f}) | {z} |"
        )
    md.append("")
    md.append(
        f"Local density (neighbours within {radius} {units}): perturbed {rep['density_by_group']['perturbed']}, NTC {rep['density_by_group']['ntc']}, unassigned {rep['density_by_group']['unassigned']:.2f}."
    )
    md.append("")
    md.append("![overview](overview.png)")
    (out / "audit.md").write_text("\n".join(md) + "\n")
    print(
        json.dumps(
            {k: rep[k] for k in ("n_units", "guide_call_rate", "n_perturbed", "n_ntc", "n_targets")}
        ),
        flush=True,
    )
    return rep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", nargs="+", required=True)
    ap.add_argument("--out", default="reports/audit")
    ap.add_argument("--n-perm", type=int, default=100)
    a = ap.parse_args()
    for d in a.dataset:
        audit(d, Path(a.out), a.n_perm)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
