"""Loader for the Perturb-FISH in vivo tumour xenograft (Binan et al., Cell 2025).

Data: Brain Image Library doi:10.35077/ace-gem-get, folder
``extras/tumors/processed/finaltables``. One FFPE section of a humanised NSG mouse tumour
xenograft imaged on a MERSCOPE (500-gene immuno-oncology panel) after Perturb-FISH readout
of a 35-target TLR/NF-kB CRISPR library (2 guides per target plus control guides).

Table linkage (verified in ``reports/audit/perturb_fish_inspection.md``):

- ``merfishcounttable.csv`` (187,215 x 553, no header): column 1 is the cell index 1..N,
  column 2 the segmentation area in mosaic pixels, column 3 the MATLAB linear index of the
  centroid in a 169,973 x 102,389 mosaic, columns 4..503 the 500 gene counts in codebook
  order, columns 504..553 fifty blank barcodes.
- ``coordinates.csv`` (187,215 x 2, no header): centroid (x = mosaic column, y = mosaic row)
  in mosaic pixels; row i belongs to count-table row i (``idx == y + (x - 1) * 169973``).
- ``tumorMerfish.csv`` header gives the 500 gene names in count-table column order.
- ``allcellsPerturbationTable.csv`` (187,215 x 77 binary, no header): guide detections per
  cell. Columns 2i-1 and 2i are the two guides of target i (i = 1..35, order below), columns
  71..74 are control guides, 75..77 are unused barcodes. The authors zero columns 7, 30 and
  68 as noisy decodes (``tumorpreprocessing/preparetablesforPython.m`` in
  github.com/lbinan/Perturb-FISH); this loader does the same and keeps them under ``obs``.

Coordinates are converted to micrometres with ``um_per_px`` (default 0.108, the MERSCOPE
mosaic pixel; the span of the author coordinates divided by the span of Vizgen's
``cell_metadata.csv`` centres gives 9.18 to 9.47 px/um, i.e. 0.106 to 0.109 um/px).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from spatialspill.schema import (
    NONE_LABEL,
    NTC_LABEL,
    SpatialScreenSchema,
    edge_distance_from_coords,
)

MOSAIC_ROWS = 169_973
MOSAIC_COLS = 102_389
N_GENES = 500
N_GUIDE_COLS = 77
NOISY_GUIDE_COLS: tuple[int, ...] = (7, 30, 68)  # 1-based, zeroed by the authors
CONTROL_GUIDE_COLS: tuple[int, ...] = (71, 72, 73, 74)  # 1-based
UNUSED_GUIDE_COLS: tuple[int, ...] = (75, 76, 77)  # 1-based

# Pooled target order used by the authors (tumorpooledperturbations.csv row names and the
# perturbation index map in tumoranalysis/clusterCells.ipynb).
TARGETS: tuple[str, ...] = (
    "CD14", "CHUK", "IKBKB", "IRAK1", "IRAK4", "IRF3", "IRF5", "IRF7", "JUN", "LBP",
    "LY96", "MAP2K1", "MAP2K2", "MAP2K3", "MAP2K4", "MAP2K6", "MAP2K7", "MAP3K7", "MAPK14",
    "MYD88", "NFKB1", "NFKBIA", "PELI1", "PIK3CA", "RELA", "RIPK1", "TAB1", "TAB2", "TBK1",
    "TICAM1", "TIRAP", "TLR4", "TRADD", "TRAF6", "TRAM1",
)  # fmt: skip

# Author quality filter from preparetablesforPython.m (counts and pixel area).
AUTHOR_FILTER = {"min_counts": 40, "max_counts": 800, "min_area_px": 2000, "max_area_px": 30000}


def guide_names() -> list[str]:
    """Names for the 77 columns of allcellsPerturbationTable.csv (1-based numbering)."""
    names: list[str] = []
    for t in TARGETS:
        names += [f"{t}_g1", f"{t}_g2"]
    names += [f"Control_g{k}" for k in range(1, len(CONTROL_GUIDE_COLS) + 1)]
    names += [f"unused_{c}" for c in UNUSED_GUIDE_COLS]
    assert len(names) == N_GUIDE_COLS
    return names


def _final_tables(raw_dir: Path) -> Path:
    return raw_dir / "extras" / "tumors" / "processed" / "finaltables"


def read_count_table(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (gene counts int32, blank counts int32, area px, centroid linear index)."""
    df = pd.read_csv(path, header=None, dtype=np.float64, engine="c")
    arr = df.to_numpy()
    if arr.shape[1] != 3 + N_GENES + 50:
        raise ValueError(f"unexpected count-table width {arr.shape[1]}")
    if not np.array_equal(arr[:, 0], np.arange(1, len(arr) + 1)):
        raise ValueError("count-table column 1 is not 1..N")
    genes = arr[:, 3 : 3 + N_GENES].astype(np.int32)
    blanks = arr[:, 3 + N_GENES :].astype(np.int32)
    return genes, blanks, arr[:, 1], arr[:, 2]


def read_gene_names(path: Path) -> list[str]:
    names = list(pd.read_csv(path, nrows=0).columns)
    if len(names) != N_GENES:
        raise ValueError(f"expected {N_GENES} gene names in {path.name}, got {len(names)}")
    return names


def _match_rows(pool: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Index in ``pool`` of each row of ``query`` by exact match, -1 if absent or ambiguous."""
    lut: dict[bytes, int] = {}
    dup: set[bytes] = set()
    for i, row in enumerate(np.ascontiguousarray(pool)):
        k = row.tobytes()
        if k in lut:
            dup.add(k)
        else:
            lut[k] = i
    out = np.full(len(query), -1, dtype=np.int64)
    for j, row in enumerate(np.ascontiguousarray(query)):
        k = row.tobytes()
        if k in lut and k not in dup:
            out[j] = lut[k]
    return out


def assign_guides(P: np.ndarray) -> pd.DataFrame:
    """Per-cell guide call from the 77-column binary table (after zeroing noisy columns)."""
    names = guide_names()
    P = P.astype(np.int16).copy()
    P[:, [c - 1 for c in NOISY_GUIDE_COLS]] = 0
    P[:, [c - 1 for c in UNUSED_GUIDE_COLS]] = 0
    n_guides = P.sum(axis=1)
    first = P.argmax(axis=1)
    guide = np.where(n_guides == 0, NONE_LABEL, np.array(names)[first])
    guide = np.where(n_guides > 1, "multi", guide)
    target = np.array(
        [g.rsplit("_g", 1)[0] if g not in (NONE_LABEL, "multi") else NONE_LABEL for g in guide]
    )
    is_ntc = target == "Control"
    target = np.where(is_ntc, NTC_LABEL, target)
    is_perturbed = (n_guides == 1) & ~is_ntc
    return pd.DataFrame(
        {
            "guide": guide,
            "target": target,
            "is_ntc": is_ntc,
            "is_perturbed": is_perturbed,
            "n_guides": n_guides.astype(int),
        }
    )


def load_perturb_fish_tumor(
    raw_dir: str | Path = "data/raw/perturb_fish",
    um_per_px: float = 0.108,
    keep_blank_probes: bool = False,
) -> AnnData:
    """Perturb-FISH tumour xenograft as a spatialspill AnnData (all 187,215 segmented cells).

    No cell filtering is applied; the authors' filter is recorded in ``obs["author_qc_pass"]``.
    """
    raw_dir = Path(raw_dir)
    ft = _final_tables(raw_dir)
    genes, blanks, area_px, lin_idx = read_count_table(ft / "merfishcounttable.csv")
    gene_names = read_gene_names(ft / "tumorMerfish.csv")
    xy_px = pd.read_csv(ft / "coordinates.csv", header=None, dtype=np.float64).to_numpy()
    n = len(genes)
    if xy_px.shape != (n, 2):
        raise ValueError(f"coordinates.csv has shape {xy_px.shape}, expected ({n}, 2)")
    resid = np.abs(lin_idx - (xy_px[:, 1] + (xy_px[:, 0] - 1) * MOSAIC_ROWS))
    if resid.max() > 1e-3:
        raise ValueError("coordinates.csv does not match the count-table centroid index")
    P = pd.read_csv(ft / "allcellsPerturbationTable.csv", header=None, dtype=np.int8).to_numpy()
    if P.shape != (n, N_GUIDE_COLS):
        raise ValueError(f"allcellsPerturbationTable.csv has shape {P.shape}")

    obs = assign_guides(P)
    obs.index = [f"cell_{i}" for i in range(1, n + 1)]
    obs["guide_confidence"] = np.nan
    obs["area"] = area_px * um_per_px**2
    obs["area_px"] = area_px
    obs["sample"] = "tumor_sample1"
    obs["batch"] = "tumor_sample1"
    xy = xy_px * um_per_px
    obs["edge_distance"] = edge_distance_from_coords(xy, obs["sample"].to_numpy())
    obs["total_counts"] = genes.sum(axis=1)
    obs["blank_counts"] = blanks.sum(axis=1)
    f = AUTHOR_FILTER
    obs["author_qc_pass"] = (
        (obs["total_counts"] > f["min_counts"])
        & (obs["total_counts"] < f["max_counts"])
        & (area_px > f["min_area_px"])
        & (area_px < f["max_area_px"])
    )

    # Author labels recovered by exact expression match (see module docstring).
    obs["cell_type"] = "unknown"
    obs["target_published"] = NONE_LABEL
    obs["has_immune_neighbor"] = "unknown"
    tcell = ft / "perturbed_tcells_merfish.csv"
    if tcell.exists():
        tc = pd.read_csv(tcell, dtype=np.int32)
        cols = [gene_names.index(g) for g in tc.columns]
        hit = _match_rows(genes[:, cols], tc.to_numpy())
        obs.loc[obs.index[hit[hit >= 0]], "cell_type"] = "T cell"
    pub = ft / "tumorMerfish.csv"
    design = ft / "tumorpooledperturbations.csv"
    if pub.exists() and design.exists():
        pm = pd.read_csv(pub, dtype=np.int32).to_numpy()
        hit = _match_rows(genes, pm)
        D = pd.read_csv(design, index_col=0).to_numpy().T  # cells x 36
        tnames = list(pd.read_csv(design, index_col=0).index)
        nz = (D > 0).sum(axis=1)
        lab = np.where(nz == 1, np.array(tnames)[D.argmax(axis=1)], "multi")
        lab = np.where(lab == "Control", NTC_LABEL, lab)
        ok = hit >= 0
        obs.loc[obs.index[hit[ok]], "target_published"] = lab[ok]
    for fname, val in (
        ("withimmuneneighborMerfish.csv", "yes"),
        ("withoutimmuneneighborMerfish.csv", "no"),
    ):
        p = ft / fname
        if p.exists():
            hit = _match_rows(genes, pd.read_csv(p, dtype=np.int32).to_numpy())
            obs.loc[obs.index[hit[hit >= 0]], "has_immune_neighbor"] = val

    X = genes.astype(np.float32)
    var_names = list(gene_names)
    if keep_blank_probes:
        X = np.hstack([X, blanks.astype(np.float32)])
        var_names += [f"Blank-{k}" for k in range(1, blanks.shape[1] + 1)]
    adata = AnnData(X=sp.csr_matrix(X), obs=obs)
    adata.var_names = var_names
    adata.var["is_blank"] = [v.startswith("Blank-") for v in var_names]
    adata.obsm["spatial"] = xy
    adata.obsm["spatial_px"] = xy_px
    adata.uns["spatialspill"] = SpatialScreenSchema(
        technology="Perturb-FISH (MERFISH)",
        dataset="BIL ace-gem-get tumors",
        organism="Homo sapiens (tumour) in humanised NSG mouse",
        units="um",
        notes={
            "source": "doi:10.35077/ace-gem-get extras/tumors/processed/finaltables",
            "um_per_px": str(um_per_px),
            "coordinate_frame": "author MERFISH mosaic (169973 x 102389 px), x = column, y = row",
            "guide_table": "allcellsPerturbationTable.csv; cols 2i-1,2i = target i; 71-74 control;"
            " 7,30,68 zeroed (noisy); 75-77 unused",
            "area": "author segmentation area in mosaic px times um_per_px^2",
            "blank_probes": "50 blank barcodes dropped" if not keep_blank_probes else "kept",
            "author_qc": "total_counts in (40, 800) and area_px in (2000, 30000)",
            "cell_type": "T cell = rows of perturbed_tcells_merfish.csv; others unknown",
            "target_published": "authors' pooled design for 9,432 perturbed cells"
            " (tumorpooledperturbations.csv)",
        },
    ).to_uns()
    return adata
