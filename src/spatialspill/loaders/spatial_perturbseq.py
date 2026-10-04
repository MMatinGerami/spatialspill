"""Loader for the in vivo Spatial Perturb-seq mouse hippocampus dataset (Shen et al., 2026).

Data: GEO GSE274447 (``GSE274447_RAW.tar``), three Stereo-seq "adjusted cellbin" GEF files
(HDF5, SAW v7, GEF version 3), one chip per mouse. An AAV sgRNA library targeting 17 genes
plus the mSafe1 safe-harbour control was injected stereotaxically; the sgRNA transcripts are
captured like genes and appear in the gene table as 18 ``sgrna_<target>`` features
(``reports/audit/spatial_perturbseq_inspection.md``).

GEF layout used here (group ``cellBin``):

- ``cell``: structured rows ``id, x, y, offset, geneCount, expCount, dnbCount, area, ...``;
  ``offset`` indexes into ``cellExp`` so ``cell.offset`` is a CSR ``indptr``.
- ``cellExp``: ``geneID, count`` triplets, cell-major.
- ``gene``: ``geneName`` (bytes), ``offset``, ``cellCount``, ``expCount``. Chips B03018A2 and
  A03599E2 contain hundreds of empty-named zero-count rows; they are dropped.

Coordinates and area are in DNB units; the root attribute ``resolution`` gives nm per DNB
(500), so 0.5 um per DNB and 0.25 um2 per DNB pixel.

Guide calls follow the authors' R package (github.com/kimberle9/spatialperturbseq,
``annotate_gRNA``): one positive sgRNA feature gives the guide, several give "multi"
(the authors write "Multiple"), none gives "none". The authors use count > 0; the
``min_guide_umi`` argument makes the threshold explicit. ``guide_confidence`` is the fraction
of a cell's guide UMIs on its dominant guide (NaN without guide UMIs).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import h5py
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

GUIDE_PREFIX = "sgrna_"
NTC_GUIDE = "sgrna_msafe"

# GEF file per chip (GSM accession prefix as deposited in GSE274447_RAW.tar).
CHIP_FILES: dict[str, str] = {
    "C02943C3": "GSM8449354_C02943C3.adjusted.cellbin.gef",  # Mouse 1 (2 sections)
    "B03018A2": "GSM9659883_B03018A2.adjusted.cellbin.gef",  # Mouse 2
    "A03599E2": "GSM9659884_A03599E2.adjusted.cellbin.gef",  # Mouse 3
}

# Authors' sgRNA feature suffix -> standard mouse gene symbol.
GUIDE_TARGETS: dict[str, str] = {
    "C9orf72": "C9orf72",
    "Cfap410": "Cfap410",
    "clu": "Clu",
    "dpp5": "Dpp6",
    "fasn": "Fasn",
    "flcn": "Flcn",
    "gfap": "Gfap",
    "lrrk2": "Lrrk2",
    "msafe": NTC_LABEL,
    "ndufaf": "Ndufaf2",
    "oligo2": "Olig2",
    "rbfox": "Rbfox3",
    "rraga": "Rraga",
    "sh3gl2": "Sh3gl2",
    "srf": "Srf",
    "stk39": "Stk39",
    "tbk1": "Tbk1",
    "trem2": "Trem2",
}


def guide_to_target(guide: str) -> str:
    """Map an ``sgrna_<suffix>`` feature name to its target symbol (``NTC`` for mSafe)."""
    suffix = guide.removeprefix(GUIDE_PREFIX)
    if suffix not in GUIDE_TARGETS:
        raise KeyError(f"unknown sgRNA feature {guide!r}")
    return GUIDE_TARGETS[suffix]


def read_cellbin_gef(path: str | Path) -> tuple[pd.DataFrame, sp.csr_matrix, np.ndarray]:
    """Read a cellbin GEF into (cells table, cells x genes CSR counts, gene names).

    ``cells`` has columns ``cell_id, x_um, y_um, area_um2, n_genes, n_umi``. Empty-named gene
    rows are dropped from the matrix and from ``genes``.
    """
    with h5py.File(path, "r") as h:
        res_nm = float(np.asarray(h.attrs["resolution"]).ravel()[0])
        cell = h["cellBin/cell"][:]
        gene = h["cellBin/gene"][:]
        exp = h["cellBin/cellExp"][:]
    um_per_dnb = res_nm / 1000.0
    genes = np.array([g.decode() for g in gene["geneName"]], dtype=object)
    n_cells = len(cell)
    indptr = np.append(cell["offset"], cell["offset"][-1] + cell["geneCount"][-1]).astype(np.int64)
    if indptr[-1] != len(exp):
        raise ValueError(f"{Path(path).name}: cell offsets do not tile cellExp")
    X = sp.csr_matrix(
        (exp["count"].astype(np.float32), exp["geneID"].astype(np.int64), indptr),
        shape=(n_cells, len(gene)),
    )
    keep = genes != ""
    if not keep.all():
        X = X[:, np.flatnonzero(keep)]
        genes = genes[keep]
    cells = pd.DataFrame(
        {
            "cell_id": cell["id"].astype(np.int64),
            "x_um": cell["x"].astype(float) * um_per_dnb,
            "y_um": cell["y"].astype(float) * um_per_dnb,
            "area_um2": cell["area"].astype(float) * um_per_dnb**2,
            "n_genes": cell["geneCount"].astype(np.int64),
            "n_umi": cell["expCount"].astype(np.int64),
        }
    )
    return cells, X, genes


def assign_guides(G: np.ndarray, guide_names: list[str], min_guide_umi: int = 1) -> pd.DataFrame:
    """Per-cell guide call from a dense cells x guides UMI matrix."""
    names = np.array(guide_names, dtype=object)
    pos = G >= min_guide_umi
    n_guides = pos.sum(axis=1)
    total = G.sum(axis=1)
    dominant = G.argmax(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        conf = np.where(total > 0, G.max(axis=1) / np.where(total > 0, total, 1), np.nan)
    guide = np.where(n_guides == 1, names[dominant], np.where(n_guides > 1, "multi", NONE_LABEL))
    target = np.array(
        [guide_to_target(g) if g.startswith(GUIDE_PREFIX) else NONE_LABEL for g in guide],
        dtype=object,
    )
    is_ntc = target == NTC_LABEL
    is_perturbed = (n_guides == 1) & ~is_ntc
    return pd.DataFrame(
        {
            "guide": guide.astype(str),
            "guide_confidence": conf.astype(float),
            "target": target.astype(str),
            "is_ntc": is_ntc,
            "is_perturbed": is_perturbed,
            "n_guides": n_guides.astype(int),
            "guide_umis_total": total.astype(int),
        }
    )


def _load_chip(path: Path, chip: str, min_guide_umi: int) -> AnnData:
    cells, X, genes = read_cellbin_gef(path)
    is_guide = np.array([g.startswith(GUIDE_PREFIX) for g in genes])
    guide_names = [str(g) for g in genes[is_guide]]
    G = X[:, np.flatnonzero(is_guide)].toarray()
    obs = assign_guides(G, guide_names, min_guide_umi)
    obs.index = pd.Index([f"{chip}_{i}" for i in cells["cell_id"]])
    obs["area"] = cells["area_um2"].to_numpy()
    obs["cell_type"] = "unknown"
    obs["sample"] = chip
    obs["batch"] = chip
    obs["n_umi"] = cells["n_umi"].to_numpy()
    obs["n_genes"] = cells["n_genes"].to_numpy()
    xy = cells[["x_um", "y_um"]].to_numpy()
    obs["edge_distance"] = edge_distance_from_coords(xy, obs["sample"].to_numpy())
    adata = AnnData(X=X[:, np.flatnonzero(~is_guide)], obs=obs)
    adata.var_names = pd.Index([str(g) for g in genes[~is_guide]])
    adata.obsm["spatial"] = xy
    adata.obsm["guide_counts"] = pd.DataFrame(
        G.astype(np.int32), index=obs.index, columns=guide_names
    )
    return adata


def load_spatial_perturbseq(
    raw_dir: str | Path = "data/raw/spatial_perturbseq",
    chips: list[str] | tuple[str, ...] | None = None,
    min_guide_umi: int = 1,
) -> AnnData:
    """Spatial Perturb-seq chips as one spatialspill AnnData (all segmented cells, no filtering).

    ``X`` holds raw counts of endogenous genes (union across chips, missing genes filled with
    0); the 18 sgRNA features live in ``obsm["guide_counts"]`` (DataFrame, one column per
    guide). GEF files are looked up in ``raw_dir/extracted`` then ``raw_dir``.
    """
    raw_dir = Path(raw_dir)
    chips = list(CHIP_FILES) if chips is None else list(chips)
    unknown = [c for c in chips if c not in CHIP_FILES]
    if unknown:
        raise ValueError(f"unknown chips {unknown}; known: {list(CHIP_FILES)}")
    parts: list[AnnData] = []
    for chip in chips:
        fname = CHIP_FILES[chip]
        path = raw_dir / "extracted" / fname
        if not path.exists():
            path = raw_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"{fname} not found under {raw_dir}")
        parts.append(_load_chip(path, chip, min_guide_umi))
    adata = parts[0] if len(parts) == 1 else ad.concat(parts, join="outer", fill_value=0)
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata.obs["batch"] = adata.obs["batch"].astype(str)
    guide_names = list(adata.obsm["guide_counts"].columns)
    adata.uns["spatialspill"] = SpatialScreenSchema(
        technology="Spatial Perturb-seq (Stereo-seq)",
        dataset="GSE274447",
        organism="Mus musculus",
        units="um",
        notes={
            "source": "GEO GSE274447_RAW.tar, SAW v7.0.0 adjusted cellbin GEF, mm10",
            "chips": ",".join(chips),
            "coordinates": "GEF cell centroid x,y in DNB units times 0.5 um (resolution attr 500 nm)",
            "area": "GEF cell area in DNB pixels times 0.25 um2",
            "guide_features": ",".join(guide_names),
            "guide_call": f"authors' annotate_gRNA rule with count >= {min_guide_umi}: "
            "one positive sgrna_* feature -> guide, several -> multi, none -> none",
            "guide_confidence": "fraction of the cell's guide UMIs on the dominant guide",
            "ntc": f"{NTC_GUIDE} (mSafe1 safe-harbour control)",
            "cell_type": "not provided by the authors (GEF cellTypeID all 0)",
            "x_excludes": "sgrna_* features (in obsm['guide_counts']) and empty-named gene rows",
        },
    ).to_uns()
    return adata
