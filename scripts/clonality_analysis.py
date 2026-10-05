"""Exploratory: why are same-target cells spatially clustered?

For each dataset, among pairs of guide-assigned cells within D_max:
  (a) distance profile of same-guide pairs, same-target-different-guide pairs, different-target
      pairs, relative to the stratified permutation null;
  (b) connected-component ("clone") size distribution of same-guide cells on the pruned
      Delaunay graph;
  (c) whether UNASSIGNED cells adjacent to cells with target g are shifted towards g's
      autonomous expression profile (misassignment / bleed-through signature): for each g with
      enough cells, correlate the autonomous profile (mean of g cells minus mean of NTC cells)
      with the neighbour profile (mean of unassigned cells adjacent to g cells minus mean of
      unassigned cells adjacent to NTC cells), both on the log-normalised scale.

Usage: uv run python scripts/clonality_analysis.py --dataset perturb_fish perturb_multi
Writes reports/audit/<dataset>/clonality.{json,md,png}.
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
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.graphs import GraphConfig, build_graph, pairwise_distances_within
from spatialspill.loaders.registry import load_dataset
from spatialspill.outcomes import lognorm
from spatialspill.permutation import permute_within_strata, strata_codes

KW: dict[str, dict] = {
    "perturb_multi": {"batches": ["4", "5", "8", "9", "10"], "max_cells_per_batch": 150_000}
}
DMAX = {"default": 100.0, "perturb_map": 400.0, "perturb_dbit": 3.5}
BINS = {
    "default": [0, 15, 30, 50, 75, 100],
    "perturb_map": [0, 110, 210, 310, 400],
    "perturb_dbit": [0, 1.5, 2.5, 3.5],
}


def pair_profile(
    D: sp.csr_matrix, guide: np.ndarray, target: np.ndarray, assigned: np.ndarray, bins: list[float]
) -> pd.DataFrame:
    Dc = sp.triu(D, k=1).tocoo()
    m = assigned[Dc.row] & assigned[Dc.col]
    r, c, d = Dc.row[m], Dc.col[m], Dc.data[m]
    same_guide = guide[r] == guide[c]
    same_target = (target[r] == target[c]) & ~same_guide
    cat = np.where(
        same_guide,
        "same_guide",
        np.where(same_target, "same_target_other_guide", "different_target"),
    )
    df = pd.DataFrame({"d": d, "cat": cat})
    df["bin"] = pd.cut(df["d"], bins, include_lowest=True)
    return df.groupby(["bin", "cat"], observed=False).size().unstack(fill_value=0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", nargs="+", required=True)
    ap.add_argument("--n-perm", type=int, default=20)
    ap.add_argument("--out", default="reports/audit")
    a = ap.parse_args()
    for name in a.dataset:
        adata = load_dataset(name, **KW.get(name, {}))
        out = Path(a.out) / name
        out.mkdir(parents=True, exist_ok=True)
        obs = adata.obs
        dmax = DMAX.get(name, DMAX["default"])
        bins = BINS.get(name, BINS["default"])
        guide = obs["guide"].astype(str).to_numpy()
        target = obs["target"].astype(str).to_numpy()
        assigned = (obs["is_perturbed"] | obs["is_ntc"]).to_numpy()
        D = pairwise_distances_within(adata, dmax)
        obs_prof = pair_profile(D, guide, target, assigned, bins)
        strata = strata_codes(obs, ("sample", "cell_type"))
        rng = np.random.default_rng(0)
        nulls = []
        for _ in range(a.n_perm):
            perm = rng.permutation(
                len(guide)
            )  # joint permutation of (guide, target) pairs within strata
            pg = permute_within_strata(np.arange(len(guide)), strata, rng)
            del perm
            nulls.append(pair_profile(D, guide[pg], target[pg], assigned[pg], bins))
        null_mean = sum(nulls) / len(nulls)
        ratio = (obs_prof / null_mean.replace(0, np.nan)).round(2)

        # clone sizes on pruned Delaunay among same-guide cells
        build_graph(
            adata,
            GraphConfig(kind="delaunay", max_edge_um=50.0 if name != "perturb_map" else 120.0),
        )
        A = adata.obsp["spatial_connectivities"].tocoo()
        keep = assigned[A.row] & assigned[A.col] & (guide[A.row] == guide[A.col])
        G = sp.coo_matrix((np.ones(keep.sum()), (A.row[keep], A.col[keep])), shape=A.shape).tocsr()
        _ncomp, lab = connected_components(G, directed=False)
        sizes = pd.Series(lab[assigned]).value_counts()
        clone_sizes = sizes.value_counts().sort_index()
        frac_in_clones = float((sizes[sizes >= 2]).sum() / assigned.sum())

        # misassignment signature: unassigned neighbours of g-cells vs of NTC cells
        sig: dict[str, float] = {}
        if (
            adata.n_vars > 500
        ):  # whole-transcriptome assays: restrict to the 500 most detected genes
            tot = np.asarray(adata.X.sum(axis=0)).ravel()
            adata = adata[:, np.sort(np.argsort(-tot)[:500])].copy()
        Y = lognorm(adata)
        A1 = adata.obsp["spatial_connectivities"].tocsr()
        unassigned = ~assigned
        ntc = obs["is_ntc"].to_numpy()
        if ntc.sum() >= 20:
            ntc_mean = Y[ntc].mean(0)
            nb_ntc = np.asarray(A1[:, ntc].sum(1)).ravel() > 0
            nb_ntc_mean = Y[unassigned & nb_ntc].mean(0)
            cors = []
            for g in pd.Series(target[obs["is_perturbed"].to_numpy()]).value_counts().index[:60]:
                gm = target == g
                if gm.sum() < 30:
                    continue
                auto = Y[gm].mean(0) - ntc_mean
                nb_g = (np.asarray(A1[:, gm].sum(1)).ravel() > 0) & unassigned & ~nb_ntc
                if nb_g.sum() < 30:
                    continue
                nb = Y[nb_g].mean(0) - nb_ntc_mean
                cors.append((g, float(np.corrcoef(auto, nb)[0, 1]), int(gm.sum()), int(nb_g.sum())))
            cdf = pd.DataFrame(
                cors, columns=["target", "cor_auto_vs_unassigned_neighbours", "n_g", "n_neighbours"]
            )
            cdf.to_csv(out / "misassignment_signature.csv", index=False)
            sig = {
                "median_cor": float(cdf["cor_auto_vs_unassigned_neighbours"].median()),
                "n_targets": len(cdf),
                "frac_cor_above_0.3": float(
                    (cdf["cor_auto_vs_unassigned_neighbours"] > 0.3).mean()
                ),
            }
        rep = {
            "dataset": name,
            "d_max": dmax,
            "n_assigned": int(assigned.sum()),
            "pair_counts_observed": {str(k): v.to_dict() for k, v in obs_prof.iterrows()},
            "pair_ratio_obs_over_null": {str(k): v.to_dict() for k, v in ratio.iterrows()},
            "same_guide_components": {
                "n_components": len(sizes),
                "fraction_assigned_in_components_ge2": frac_in_clones,
                "size_distribution": {int(k): int(v) for k, v in clone_sizes.items()},
            },
            "misassignment_signature": sig,
        }
        (out / "clonality.json").write_text(json.dumps(rep, indent=2, default=str) + "\n")
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ratio.plot(kind="bar", ax=ax[0], logy=True)
        ax[0].set_ylabel("observed / null pair count")
        ax[0].set_title(f"{name}: pairs within {dmax} {adata.uns['spatialspill']['units']}")
        clone_sizes.plot(kind="bar", ax=ax[1], logy=True, color="grey")
        ax[1].set_xlabel("same-guide component size (pruned Delaunay)")
        ax[1].set_title(f"{frac_in_clones:.1%} of assigned cells in components >= 2")
        fig.tight_layout()
        fig.savefig(out / "clonality.png", dpi=120)
        md = [
            f"# Same-target clustering: {name}",
            "",
            "Exploratory (not pre-registered). Generated by scripts/clonality_analysis.py.",
            "",
            "## Observed / null pair counts by distance bin",
            "",
            ratio.to_markdown(),
            "",
            f"## Same-guide connected components (pruned Delaunay): {frac_in_clones:.1%} of assigned cells in components of size >= 2",
            "",
            clone_sizes.to_frame("n_components").T.to_markdown(),
            "",
            f"## Misassignment signature: {json.dumps(sig)}",
            "",
            "![clonality](clonality.png)",
        ]
        (out / "clonality.md").write_text("\n".join(md) + "\n")
        print(
            name,
            json.dumps(
                {"ratio": rep["pair_ratio_obs_over_null"], "clones": frac_in_clones, "sig": sig},
                default=str,
            )[:1500],
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
