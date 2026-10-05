"""Downstream analysis of one real-data estimate table: artifact fraction (C3), niche gene
ranking and ligand-receptor enrichment (C6), held-out-target prediction from embeddings (C5).

Usage: uv run python scripts/analyze_real.py --results results/<hash> [--embeddings ...]
Writes into the same results directory: artifact_fraction.json, bleedthrough_table.csv,
niche_ranking.csv, lr_enrichment.json, prediction.json, prediction_per_target.csv.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spatialspill.artifacts import artifact_fraction, bleedthrough_table
from spatialspill.biology import load_lr_pairs, lr_enrichment, rank_niche_genes
from spatialspill.embeddings import go_embedding, perturbseq_embedding
from spatialspill.prediction import loto_predict, response_matrix

SUPPORT = Path(__file__).resolve().parents[1] / "data" / "raw" / "supporting"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument(
        "--rings",
        nargs="*",
        type=int,
        default=[1, 2],
        help="rings used for spillover profiles (ring 0 excluded by default: calibration not established)",
    )
    ap.add_argument("--organism", default="human", choices=["human", "mouse"])
    a = ap.parse_args()
    d = Path(a.results)
    files = glob.glob(str(d / "*_e2_estimates.csv")) + glob.glob(str(d / "*_e1_estimates.csv"))
    if not files:
        raise SystemExit("no estimate table found")
    df = pd.read_csv(files[0])
    out: dict[str, object] = {"table": Path(files[0]).name, "rings_used": a.rings}

    # C3: bleed-through fraction (first ring vs last ring), on the pooled cell-type group
    bt = bleedthrough_table(df)
    bt.to_csv(d / "bleedthrough_table.csv", index=False)
    af = artifact_fraction(bt)
    out["artifact_fraction"] = af
    (d / "artifact_fraction.json").write_text(json.dumps(af, indent=2) + "\n")

    # C6: ranking and LR enrichment (rings as given)
    sp = df[(df["kind"] == "spillover") & df["ring"].isin(a.rings)]
    rk = rank_niche_genes(sp, df)
    rk.to_csv(d / "niche_ranking.csv", index=False)
    lr = load_lr_pairs(
        omnipath_tsv=SUPPORT / "lr" / "omnipath_ligrec.tsv",
        cellphonedb_csv=SUPPORT / "lr" / "cellphonedb_interaction_input.csv",
    )
    enr = lr_enrichment(sp[sp["cell_type"] == "all"], lr)
    out["lr_enrichment"] = enr
    (d / "lr_enrichment.json").write_text(json.dumps(enr, indent=2) + "\n")

    # C5: held-out-target prediction of autonomous and spillover z-profiles
    pred: dict[str, object] = {}
    R_auto = response_matrix(df, "autonomous")
    R_spill = response_matrix(df, "spillover", rings=a.rings)
    R_auto.index = R_auto.index.astype(str).str.upper()
    R_spill.index = R_spill.index.astype(str).str.upper()
    embs: dict[str, pd.DataFrame] = {}
    gwps = SUPPORT / "replogle" / "K562_gwps_normalized_bulk_01.h5ad"
    if gwps.exists():
        embs["replogle_k562_gwps"] = perturbseq_embedding(gwps)
    gaf = SUPPORT / "go" / ("goa_human.gaf.gz" if a.organism == "human" else "goa_mouse.gaf.gz")
    if gaf.exists():
        embs["go"] = go_embedding(gaf, SUPPORT / "go" / "go-basic.obo")
    rows = []
    for ename, E in embs.items():
        for kind, R in (("autonomous", R_auto), ("spillover", R_spill)):
            if R.empty:
                continue
            res = loto_predict(R, E, n_null=20)
            per = res.pop("per_target")
            if isinstance(per, pd.DataFrame) and len(per):
                per = per.assign(embedding=ename, kind=kind)
                rows.append(per)
            pred[f"{ename}:{kind}"] = res
    out["prediction"] = pred
    (d / "prediction.json").write_text(json.dumps(pred, indent=2, default=float) + "\n")
    if rows:
        pd.concat(rows).to_csv(d / "prediction_per_target.csv", index=False)
    (d / "analysis_summary.json").write_text(json.dumps(out, indent=2, default=float) + "\n")
    print(json.dumps(out, indent=2, default=float))
    print(rk.head(15).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
