"""Loader for Perturb-DBiT (Baysoy et al., Nat Biotechnol 2026; GEO GSE319277).

Unit of analysis: the microfluidic pixel (not the cell). Coordinates are pixel grid indices
(barcode A x barcode B) multiplied by ``pitch_um``; the pitch is not stated in the GEO
metadata, so the default is 1.0 and distances are in pixel units unless a pitch is supplied.
A pixel's guide call is the dominant sgRNA by UMI count (see ``dominant_guide``).

Samples with both an expression matrix and an sgRNA table in GEO: Bru.LC.6 (total RNA),
mTSG.LV.2, mTSG.LV.3, Bru.LC.8, Brie.LU.1, Bru.LC.1.33.3. Expression TSVs come in both
orientations (pixels x genes or genes x pixels) and are detected by the ``AxB`` pattern.
"""

from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from spatialspill.loaders._util import dominant_guide, guide_to_target, is_ntc_guide
from spatialspill.schema import SpatialScreenSchema, edge_distance_from_coords

PIX_RE = re.compile(r"^(?:bc)?(\d+)x(\d+)$")

SAMPLES: dict[str, dict[str, str]] = {
    # sample: (expression file, sgRNA file) inside GSE319277_RAW.tar
    "Bru.LC.6": {
        "expr": "GSM9514752_Bru.LC.6_allRNA.tsv.gz",
        "sg": "GSM9514752_sgCounts_Bru.LC.6_grnalib.txt.gz",
        "tissue": "HT29 lung metastasis (NSG mouse), FFPE",
        "library": "Brunello-type (human)",
    },
    "mTSG.LV.2": {
        "expr": "GSM9514753_mTSG.LV.2.tsv.gz",
        "sg": "GSM9514753_sgCounts_mTSG.LV.2.txt.gz",
        "tissue": "AAV GEMM liver tumour, FFPE",
        "library": "mouse TSG",
    },
    "mTSG.LV.3": {
        "expr": "GSM9514754_mTSG.LV.3.updated.tsv.gz",
        "sg": "GSM9514754_sgCounts_mTSG.LV.3.txt.gz",
        "tissue": "AAV GEMM liver tumour, FFPE",
        "library": "mouse TSG",
    },
    "Bru.LC.8": {
        "expr": "GSM9514756_Bru.LC.8.tsv.gz",
        "sg": "GSM9514756_sgCounts_Bru.LC.8_grnalib.txt.gz",
        "tissue": "HT29 lung metastasis (NSG mouse), FFPE",
        "library": "Brunello-type (human)",
    },
    "Brie.LU.1": {
        "expr": "GSM9514757_Brie.LU.1.tsv.gz",
        "sg": "GSM9514757_sgCounts_Brie.LU.1_1.txt.gz",
        "tissue": "E0771 lung metastasis (syngeneic), FF",
        "library": "Brie (mouse)",
    },
    "Bru.LC.1.33.3": {
        "expr": "GSM9514758_Bru.LC.1.33.3.tsv.gz",
        "sg": "GSM9514758_sgCounts_Bru.LC.1.33.3.txt.gz",
        "tissue": "HT29 lung metastasis (NSG mouse), FFPE",
        "library": "Brunello-type (human)",
    },
}


def _extract(raw_dir: Path) -> Path:
    out = raw_dir / "extracted"
    tar = raw_dir / "GSE319277_RAW.tar"
    if not out.exists() or not any(out.iterdir()):
        out.mkdir(exist_ok=True)
        with tarfile.open(tar) as tf:
            tf.extractall(out, filter="data")
    return out


def _parse_pixel(s: str) -> tuple[int, int] | None:
    m = PIX_RE.match(str(s).strip().strip('"'))
    return (int(m.group(1)), int(m.group(2))) if m else None


def read_sg_table(path: Path) -> pd.DataFrame:
    """Return long table with columns pixel (``AxB``), grna, count from any of the GEO formats."""
    with gzip.open(path, "rt") as fh:
        first = fh.readline()
    if first.startswith("BA\tBB"):
        df = pd.read_csv(path, sep="\t")
        df["pixel"] = df["BA"].astype(int).astype(str) + "x" + df["BB"].astype(int).astype(str)
        return df[["pixel", "grna", "count"]]
    if first.startswith('"","barcode"'):
        df = pd.read_csv(path)
        df["pixel"] = df["barcode"].astype(str)
        return df[["pixel", "grna", "count"]]
    df = pd.read_csv(path, sep="\t", header=None)
    df = df.iloc[:, :3]
    df.columns = ["pixel", "grna", "count"]
    df["pixel"] = df["pixel"].astype(str).str.replace("^bc", "", regex=True)
    return df


def read_expression(path: Path) -> pd.DataFrame:
    """Pixels x genes count table, oriented by detecting ``AxB`` pixel ids."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    idx_is_pixel = all(_parse_pixel(i) is not None for i in list(df.index[:20]))
    if not idx_is_pixel:
        df = df.T
    df.index = [i.strip().strip('"') for i in df.index.astype(str)]
    df.columns = [re.sub(r"^X(?=\d)", "", str(c)) for c in df.columns]
    return df


def load_perturb_dbit(
    raw_dir: str | Path = "data/raw/perturb_dbit",
    samples: list[str] | None = None,
    pitch_um: float = 1.0,
    min_guide_umis: int = 1,
    min_guide_frac: float = 0.5,
    max_genes: int | None = None,
) -> AnnData:
    raw_dir = Path(raw_dir)
    ext = _extract(raw_dir)
    samples = samples or list(SAMPLES)
    parts: list[AnnData] = []
    for s in samples:
        info = SAMPLES[s]
        expr = read_expression(ext / info["expr"])
        sg = read_sg_table(ext / info["sg"])
        calls = dominant_guide(sg, "pixel", "grna", "count", min_guide_umis, min_guide_frac)
        obs = pd.DataFrame(index=expr.index)
        obs = obs.join(calls, how="left")
        obs["guide"] = obs["guide"].fillna("none").astype(str)
        for c in ("guide_umis", "guide_umis_total", "n_guides"):
            obs[c] = obs[c].fillna(0).astype(int)
        obs["target"] = [guide_to_target(g) if g != "none" else "none" for g in obs["guide"]]
        obs["is_ntc"] = np.array(
            [g != "none" and is_ntc_guide(g) for g in obs["guide"]], dtype=bool
        )
        obs["is_perturbed"] = (obs["guide"] != "none").to_numpy() & ~obs["is_ntc"].to_numpy()
        ab = np.array([_parse_pixel(i) for i in obs.index], dtype=float)
        xy = ab * pitch_um
        obs["area"] = np.nan
        obs["cell_type"] = "pixel"
        obs["sample"] = s
        obs["batch"] = s
        obs["edge_distance"] = edge_distance_from_coords(xy, obs["sample"].to_numpy())
        obs["pixel_a"] = ab[:, 0].astype(int)
        obs["pixel_b"] = ab[:, 1].astype(int)
        obs["total_counts"] = expr.sum(axis=1).to_numpy()
        X = sp.csr_matrix(expr.to_numpy(dtype=np.float32))
        ad = AnnData(X=X, obs=obs)
        ad.var_names = list(expr.columns)
        ad.var_names_make_unique()
        ad.obs_names = [f"{s}:{i}" for i in obs.index]
        ad.obsm["spatial"] = xy
        parts.append(ad)
    if len(parts) == 1:
        adata = parts[0]
    else:
        import anndata as an

        adata = an.concat(parts, join="outer", fill_value=0, merge="same")
    if max_genes is not None and adata.n_vars > max_genes:
        tot = np.asarray(adata.X.sum(axis=0)).ravel()
        keep = np.argsort(-tot)[:max_genes]
        adata = adata[:, np.sort(keep)].copy()
    adata.uns["spatialspill"] = SpatialScreenSchema(
        technology="Perturb-DBiT",
        dataset="GSE319277",
        organism="Mus musculus (host) / Homo sapiens (HT29 tumour)",
        units="pixel" if pitch_um == 1.0 else "um",
        notes={
            "unit": "microfluidic pixel",
            "pitch_um": str(pitch_um),
            "samples": ",".join(samples),
        },
    ).to_uns()
    return adata
