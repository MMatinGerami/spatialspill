"""Loader for Perturb-map Visium data (Dhainaut et al., Cell 2022; GEO GSE193460).

Unit of analysis: the Visium spot (55 um capture area, 100 um centre-to-centre on a hexagonal
array), not the cell. Four sections of KP (Kras G12D; Trp53 -/-) lung tumours carrying
Pro-Code barcoded CRISPR knockouts were profiled with 10x Visium; the GEO submission provides
Space Ranger outputs plus a per-spot annotation table.

Perturbation labels
-------------------
Spots were not genotyped individually. The authors clustered tumour spots (``leiden_clusters``)
and assigned each lesion cluster a Pro-Code phenotype by matching to imaging of the adjacent
section. The ``phenotypes`` column therefore contains lesion-level labels of three kinds:

- ``"<Gene>_<section>"`` or ``"<Gene>_<section>-<k>"`` (``Tgfbr2_1``, ``Jak2_1``, ``Ifngr2_2``,
  ``Tgfbr2_4-1``): lesion with a named knockout. ``target`` is the gene symbol and
  ``is_perturbed`` is True.
- ``"KP_<section>-<k>"``: tumour lesion without an assigned knockout. ``target`` is
  ``"KP_unlabeled"``; ``is_perturbed`` and ``is_ntc`` are both False because the knockout is
  unknown, not absent.
- ``"periphery"``: tumour periphery cluster; ``target`` is ``"none"``.

Spots with no label (``NA``) get ``guide == "none"``. There are no non-targeting control
labels in this dataset, so ``is_ntc`` is False everywhere and no spot carries target ``"NTC"``.
``guide_confidence`` and ``area`` are NaN (not provided).

Coordinates
-----------
``tissue_positions_list.csv.gz`` (Space Ranger v1 layout, no header) gives
``pxl_row_in_fullres`` and ``pxl_col_in_fullres``. The GEO images are downsampled, so the
pixel scale must be recovered from the known geometry. Two anchors are available:

- ``"pitch"`` (default): the median nearest-neighbour distance between spot centres in pixel
  space (over all 4992 array positions) equals 100 um on the Visium array. In all four
  sections this is 49.0 px, so 1 um = 0.490 px.
- ``"spot_diameter"``: ``spot_diameter_fullres`` from ``scalefactors_json`` taken as 55 um
  (32.3 px, so 1 um = 0.588 px). In Space Ranger output the ratio spot diameter / pitch is
  0.656 rather than 0.55, meaning the reported diameter corresponds to about 65 um. Using it
  with 55 um shrinks the array pitch to about 84 um, so this anchor is kept only for
  comparison.

``obsm["spatial"]`` is ``(pxl_col, pxl_row) / px_per_um``, i.e. x = image column, y = image
row, in micrometres. Both scale estimates are recorded in ``uns["spatialspill"]["notes"]``.

Spot filtering: only ``in_tissue == 1`` spots (the filtered matrix) are used, and spots listed
in ``*_off_tissue.csv.gz`` (authors' manual removal) are dropped. After this the barcodes equal
the rows of the annotation table exactly.
"""

from __future__ import annotations

import gzip
import json
import re
import warnings
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.spatial import cKDTree

from spatialspill.schema import NONE_LABEL, SpatialScreenSchema, edge_distance_from_coords

VISIUM_PITCH_UM = 100.0
VISIUM_SPOT_DIAMETER_UM = 55.0
KP_UNLABELED = "KP_unlabeled"
PERIPHERY = "periphery"

PHENOTYPE_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*)_(\d+)(?:-(\d+))?$")

SAMPLES: dict[str, dict[str, str]] = {
    # section: GEO sample accession and mouse id (Dhainaut et al. 2022, GSE193460)
    "KP_1": {"gsm": "GSM5808054", "mouse": "4.1"},
    "KP_2": {"gsm": "GSM5808055", "mouse": "4.3"},
    "KP_3": {"gsm": "GSM5808056", "mouse": "5.1"},
    "KP_4": {"gsm": "GSM5808057", "mouse": "7.2"},
}

POSITIONS_COLUMNS = ("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")


def phenotype_to_target(phenotype: str) -> str:
    """Target gene symbol from a Perturb-map phenotype label.

    ``Tgfbr2_1`` and ``Tgfbr2_4-1`` give ``Tgfbr2``; ``KP_1-1`` gives ``KP_unlabeled``;
    ``periphery`` and ``none`` give ``none``.
    """
    s = str(phenotype).strip()
    if s in (NONE_LABEL, PERIPHERY, ""):
        return NONE_LABEL
    m = PHENOTYPE_RE.match(s)
    if m is None:
        return s
    gene = m.group(1)
    return KP_UNLABELED if gene == "KP" else gene


def read_positions(path: Path) -> pd.DataFrame:
    """Space Ranger ``tissue_positions_list.csv`` (v1, headerless) indexed by barcode."""
    with gzip.open(path, "rt") as fh:
        first = fh.readline()
    header = 0 if first.startswith("barcode") else None
    df = pd.read_csv(path, header=header)
    df.columns = list(POSITIONS_COLUMNS)
    return df.set_index("barcode")


def read_scalefactors(path: Path) -> dict[str, float]:
    with gzip.open(path, "rt") as fh:
        return json.load(fh)


def pitch_px(xy_px: np.ndarray) -> float:
    """Median nearest-neighbour distance between spot centres in pixels."""
    d, _ = cKDTree(xy_px).query(xy_px, k=2)
    return float(np.median(d[:, 1]))


def _load_section(
    ext: Path, section: str, micron_scale: Literal["pitch", "spot_diameter"]
) -> tuple[AnnData, dict[str, str]]:
    import scanpy as sc

    info = SAMPLES[section]
    prefix = f"{info['gsm']}_{section}_"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # duplicate var names are made unique below
        ad = sc.read_10x_h5(ext / f"{prefix}filtered_feature_bc_matrix.h5")
    ad.var_names_make_unique()
    pos = read_positions(ext / f"{prefix}tissue_positions_list.csv.gz")
    sf = read_scalefactors(ext / f"{prefix}scalefactors_json.json.gz")
    ann = pd.read_csv(ext / f"{prefix}spot_annotation.csv.gz").set_index("barcode")
    off = pd.read_csv(ext / f"{prefix}off_tissue.csv.gz")
    off_barcodes = set(off.iloc[:, 0].astype(str))

    pos = pos.loc[ad.obs_names]
    keep = (pos["in_tissue"] == 1).to_numpy() & ~pos.index.isin(off_barcodes)
    ad = ad[keep].copy()
    pos = pos.loc[ad.obs_names]
    ann = ann.reindex(ad.obs_names)

    all_xy_px = read_positions(ext / f"{prefix}tissue_positions_list.csv.gz")[
        ["pxl_col", "pxl_row"]
    ].to_numpy(dtype=float)
    px_per_um_pitch = pitch_px(all_xy_px) / VISIUM_PITCH_UM
    px_per_um_diam = float(sf["spot_diameter_fullres"]) / VISIUM_SPOT_DIAMETER_UM
    px_per_um = px_per_um_pitch if micron_scale == "pitch" else px_per_um_diam
    xy = pos[["pxl_col", "pxl_row"]].to_numpy(dtype=float) / px_per_um

    phen = ann["phenotypes"].astype("string")
    guide = phen.fillna(NONE_LABEL).astype(str).to_numpy()
    target = np.array([phenotype_to_target(g) for g in guide], dtype=object)
    is_perturbed = ~np.isin(target, [NONE_LABEL, KP_UNLABELED])
    region = np.where(
        guide == NONE_LABEL, "unlabeled", np.where(guide == PERIPHERY, PERIPHERY, "lesion")
    )
    kmeans = ann["kmeans"].fillna("unknown").astype(str).to_numpy()

    obs = pd.DataFrame(index=ad.obs_names)
    obs["guide"] = guide
    obs["guide_confidence"] = np.nan
    obs["target"] = target.astype(str)
    obs["is_ntc"] = np.zeros(len(obs), dtype=bool)
    obs["is_perturbed"] = is_perturbed.astype(bool)
    obs["area"] = np.nan
    obs["cell_type"] = np.char.add(np.char.add(kmeans.astype(str), "_"), region.astype(str))
    obs["sample"] = section
    obs["batch"] = f"mouse_{info['mouse']}"
    obs["edge_distance"] = edge_distance_from_coords(xy, obs["sample"].to_numpy())
    obs["n_counts"] = np.asarray(ad.X.sum(axis=1)).ravel().astype(int)
    obs["n_features"] = np.asarray((ad.X > 0).sum(axis=1)).ravel().astype(int)
    leiden = ann["leiden_clusters"]
    obs["leiden_clusters"] = np.where(
        leiden.isna(), NONE_LABEL, leiden.fillna(-1).astype(int).astype(str)
    )
    obs["phenotype_raw"] = phen.fillna("NA").astype(str).to_numpy()
    obs["in_tissue"] = pos["in_tissue"].to_numpy().astype(int)
    obs["array_row"] = pos["array_row"].to_numpy().astype(int)
    obs["array_col"] = pos["array_col"].to_numpy().astype(int)

    out = AnnData(X=sp.csr_matrix(ad.X, dtype=np.float32), obs=obs, var=ad.var[["gene_ids"]])
    out.obs_names = [f"{section}:{b}" for b in ad.obs_names]
    out.obsm["spatial"] = xy
    notes = {
        f"{section}_px_per_um_pitch": f"{px_per_um_pitch:.4f}",
        f"{section}_px_per_um_spot_diameter": f"{px_per_um_diam:.4f}",
        f"{section}_mouse": info["mouse"],
        f"{section}_off_tissue_dropped": str(len(off_barcodes)),
    }
    return out, notes


def load_perturb_map(
    raw_dir: str | Path = "data/raw/perturb_map",
    samples: list[str] | None = None,
    micron_scale: Literal["pitch", "spot_diameter"] = "pitch",
) -> AnnData:
    """Load the Perturb-map Visium sections into one AnnData following the spatialspill schema.

    Parameters
    ----------
    raw_dir
        Directory containing ``extracted/`` with the GSE193460 Space Ranger outputs.
    samples
        Sections to load (subset of ``KP_1`` to ``KP_4``); all four by default.
    micron_scale
        How to convert pixel coordinates to micrometres; see the module docstring.
    """
    ext = Path(raw_dir) / "extracted"
    samples = samples or list(SAMPLES)
    parts: list[AnnData] = []
    notes: dict[str, str] = {
        "unit": "Visium spot",
        "label_level": "lesion cluster (Pro-Code phenotype matched to adjacent section)",
        "ntc": "no non-targeting control labels; is_ntc is False for every spot",
        "micron_scale": micron_scale,
        "visium_pitch_um": str(VISIUM_PITCH_UM),
        "samples": ",".join(samples),
    }
    for s in samples:
        ad, n = _load_section(ext, s, micron_scale)
        parts.append(ad)
        notes.update(n)
    if len(parts) == 1:
        adata = parts[0]
    else:
        import anndata as an

        adata = an.concat(parts, join="outer", fill_value=0, merge="same")
    adata.uns["spatialspill"] = SpatialScreenSchema(
        technology="Perturb-map (Visium)",
        dataset="GSE193460",
        organism="Mus musculus",
        units="um",
        notes=notes,
    ).to_uns()
    return adata
