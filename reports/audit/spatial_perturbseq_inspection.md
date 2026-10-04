# Spatial Perturb-seq (GSE274447) inspection

Dataset: Shen et al., "Spatial perturb-seq: single-cell functional genomics within intact tissue architecture", Nat Commun 2026, doi:10.1038/s41467-026-69677-6 (PMID 41723140, from the GEO SOFT record). In vivo AAV sgRNA library (18 loci incl. the mSafe1 safe-harbour control) injected into mouse hippocampus, Stereo-seq (SAW v7.0.0, mm10), cell segmentation ("adjusted.cellbin", Deepcell for mice 2 and 3).

Files (tar 2142842880 bytes, verified complete; extracted to `data/raw/spatial_perturbseq/extracted/`):

| GSM | Chip | Title | Size |
|---|---|---|---|
| GSM8449354 | C02943C3 | Mouse 1 (2 sections) | 253 MB |
| GSM9659883 | B03018A2 | Mouse 2 | 1.72 GB |
| GSM9659884 | A03599E2 | Mouse 3 | 169 MB |

All numbers below come from h5py/numpy code run on 2026-10-04 (scratch scripts `read_gef.py`, logs `guide_stats.log`, `extra_checks.log`).

## 1. GEF layout (HDF5, geftool v1.1.x, GEF version 3, bin_type CellBin)

Root attrs: `bin_type=CellBin`, `omics=Transcriptomics`, `resolution=500` (nm per DNB), `offsetX/offsetY` (0/3, 0/0, 0/2), `version=3`.

Group `cellBin/`:

| Dataset | dtype | Meaning |
|---|---|---|
| `cell` (n_cells,) | `id u4, x i4, y i4, offset u4, geneCount u2, expCount u2, dnbCount u2, area u2, cellTypeID u2, clusterID u2` | one row per cell; `offset` indexes into `cellExp`; `cellTypeID`/`clusterID` all 0 (`cellTypeList=['default']`), i.e. no annotation shipped |
| `cellExp` (nnz,) | `geneID u4, count u2` | cell-major triplets (CSR with `cell.offset` as indptr) |
| `gene` (n_genes,) | `geneName S64, offset u4, cellCount u4, expCount u4, maxMIDcount u2` | gene table; `offset` indexes into `geneExp` |
| `geneExp` (nnz,) | `cellID u4, count u2` | the same matrix, gene-major (CSC) |
| `cellBorder` (n_cells, 32, 2) int16 | polygon vertices relative to the cell centroid, padded with 32767 |
| `cellExon`, `cellExpExon`, `geneExon`, `geneExpExon` | exon-only counts, same shapes |
| `blockIndex`, `blockSize` | spatial block index for viewers (not needed) |

Verified for all three files: CSR row sums equal `cell.expCount`, nnz per row equals `cell.geneCount`, `cell.id` is contiguous 0..n-1. Mouse 2 and Mouse 3 have 832 and 1431 gene rows with an empty `geneName` and zero counts (SAW placeholder rows); they should be dropped when loading.

## 2. Cells, genes, coordinates

| File | Cells | Genes (named) | nnz | Median UMI/cell | Median genes/cell | Extent x (um) | Extent y (um) | Median area (um2) | Median NN centroid dist (um) |
|---|---|---|---|---|---|---|---|---|---|
| Mouse 1 C02943C3 | 39,862 | 26,280 | 15,264,103 | 596 | 317 | 855.5–10945.0 | 2254.0–12341.0 | 281.8 | 15.3 |
| Mouse 2 B03018A2 | 123,921 | 29,293 (28,461) | 106,216,307 | 1,403 | 804 | 1370.0–11401.5 | 223.0–10256.5 | 292.0 | 15.1 |
| Mouse 3 A03599E2 | 65,992 | 24,296 (22,865) | 9,760,113 | 251 | 134 | 2468.5–10130.0 | 1106.5–10264.5 | 288.0 | 15.3 |

Coordinate units: `x`, `y` and `area` are in DNB units (`resolution=500` nm, so 0.5 um per DNB and 0.25 um2 per DNB pixel). The ~10 mm extent and ~15 um nearest-neighbour spacing after conversion are consistent with a coronal hemisphere at cellular resolution. The `offsetX/offsetY` attrs are 0–3 DNBs (negligible) and relate stored coordinates to the raw chip frame. Mouse 1 is "2 sections" on one chip (two tissue pieces in one coordinate frame).

## 3. Guide (sgRNA) features

Yes. Each file contains 18 features named `sgrna_<target>` in the gene table, counted exactly like genes (the AAV sgRNA transcript is captured by the polyA Stereo-seq chemistry): `sgrna_C9orf72, sgrna_Cfap410, sgrna_clu, sgrna_dpp5, sgrna_fasn, sgrna_flcn, sgrna_gfap, sgrna_lrrk2, sgrna_msafe, sgrna_ndufaf, sgrna_oligo2, sgrna_rbfox, sgrna_rraga, sgrna_sh3gl2, sgrna_srf, sgrna_stk39, sgrna_tbk1, sgrna_trem2`. Note the authors' spellings: `dpp5` (target Dpp6), `oligo2` (Olig2), `ndufaf` (Ndufaf2), `rbfox` (Rbfox3), `msafe` (mSafe1 safe-harbour control). All 17 endogenous target genes are present under standard symbols (Trem2, Rraga, ..., Rbfox3).

How the authors assign guides (github.com/kimberle9/spatialperturbseq, `R/functions.R`): Seurat object built from the Stereopy h5ad export (feature names become `sgrna-...` because Seurat replaces `_` with `-`). `annotate_pos_neg()` labels a cell "Pos" if any sgRNA feature has count > 0, else "Neg". `annotate_gRNA()` labels each cell with the sgRNA name if exactly one sgRNA feature has count > 0, "Multiple" if more than one, "Neg" if none. There is no UMI threshold and no deconvolution: a single UMI of a guide is enough to call a cell positive. `pos_neg_DEG()` and `dotplot_upregulated()` simply drop features matching `^sgrna` from DE results.

```r
annotate_gRNA <- function(obj, gRNA_list){
  cells <- FetchData(object = obj, vars = gRNA_list)
  cells <- cells %>% mutate(GeneExpressionStatus = apply(cells, 1, function(row) {
      if (all(row == 0)) "Neg"
      else if (sum(row > 0) == 1) names(row)[which(row > 0)]
      else "Multiple" }))
  obj$gRNA <- cells$GeneExpressionStatus; obj }
```

### Guide statistics (any count > 0, same rule as the authors)

| File | Guide+ cells | % of cells | Total guide UMI | Guide UMI per positive cell | Distinct guides per positive cell | "Multiple" (authors' rule) |
|---|---|---|---|---|---|---|
| Mouse 1 | 176 / 39,862 | 0.44% | 288 | 1: 104, 2: 40, 3: 24, 4: 8 | 1: 165, 2: 9, 3: 2 | 11 |
| Mouse 2 | 4,257 / 123,921 | 3.44% | 8,575 | 1: 2413, 2: 911, 3: 394, 4: 212, 5: 121, 6: 66, 7: 41, 8: 35, 9: 25, 10+: 39 (max 24) | 1: 3203, 2: 606, 3: 251, 4: 104, 5: 49, 6: 27, 7–10: 17 | 1,054 |
| Mouse 3 | 322 / 65,992 | 0.49% | 588 | 1: 171, 2: 80, 3: 43, 4: 21, 5: 3, 6: 2, 8: 1, 9: 1 | 1: 301, 2: 17, 3: 4 | 21 |

Per-guide cells with count > 0 in Mouse 2 (cells, UMI): clu 776/1093, fasn 735/1016, stk39 665/1027, oligo2 609/901, trem2 480/762, flcn 366/468, sh3gl2 366/473, gfap 322/469, rraga 298/375, Cfap410 251/339, ndufaf 224/304, msafe 204/282, rbfox 181/226, tbk1 170/217, dpp5 159/196, srf 140/171, C9orf72 94/134, lrrk2 85/122. Mouse 1 and 3 have 2–32 and 6–67 cells per guide respectively.

Spatial pattern of guide+ cells (Mouse 2): median distance from a guide+ cell to the nearest other guide+ cell 22.5 um, 64% within 30 um (injection-site clustering). Mouse 1: 42.0 um / 35%; Mouse 3: 49.7 um / 27%.

Quick spillover-relevant check (Mouse 2, 189 "confident" cells with exactly one guide at >= 3 UMI): among neighbours within 25 um (781 cells), 27.7% carry some guide and 4.74% carry the *same* guide, versus 2.47% expected per non-matching guide. Same-guide enrichment persists at 15 um (5.07% vs 2.60%) and 40 um (4.40% vs 2.66%). 2,413 of 4,257 guide+ cells in Mouse 2 have exactly 1 guide UMI. This is the signature one would expect from a mixture of true co-infection, diffusion/segmentation spill of sgRNA transcripts, and local AAV spread; it is the quantity a spillover model must disentangle.

## 4. Loader (h5py/numpy/scipy/pandas only; validated on all three files)

```python
import h5py, numpy as np, pandas as pd, scipy.sparse as sp


def read_cellbin_gef(path):
    """Return (cells DataFrame, X csr cells x genes, genes ndarray[str])."""
    with h5py.File(path, "r") as h:
        res_nm = float(h.attrs["resolution"][0])  # 500 nm per DNB
        um_per_dnb = res_nm / 1000.0  # 0.5 um
        cell = h["cellBin/cell"][:]  # id,x,y,offset,geneCount,expCount,...
        gene = h["cellBin/gene"][:]  # geneName,offset,cellCount,expCount
        exp = h["cellBin/cellExp"][:]  # geneID,count (cell-major, via cell.offset)
    genes = np.array([g.decode() for g in gene["geneName"]])
    indptr = np.append(cell["offset"], cell["offset"][-1] + cell["geneCount"][-1]).astype(np.int64)
    assert indptr[-1] == len(exp)
    X = sp.csr_matrix(
        (exp["count"].astype(np.int32), exp["geneID"].astype(np.int64), indptr),
        shape=(len(cell), len(gene)),
    )
    cells = pd.DataFrame(
        {
            "cell_id": cell["id"].astype(np.int64),
            "x_um": cell["x"] * um_per_dnb,
            "y_um": cell["y"] * um_per_dnb,
            "area_um2": cell["area"] * um_per_dnb**2,
            "n_genes": cell["geneCount"].astype(np.int64),
            "n_umi": cell["expCount"].astype(np.int64),
        }
    )
    return cells, X, genes
    # downstream: keep = genes != ""; X = X[:, keep]; genes = genes[keep]
    # guides: gidx = np.flatnonzero(np.char.startswith(genes, "sgrna_"))
```

Load time: 0.1 s / 0.7 s / <0.1 s (Mouse 1/2/3) on this machine; Mouse 2 needs about 2 GB RAM for the triplets.

## 5. Usability for cell-level spillover analysis

Usable, with Mouse 2 as the primary chip. The dataset supplies the three required ingredients per cell: guide identity (18 `sgrna_*` features counted per segmented cell), coordinates (0.5 um resolution, converted to um, plus cell polygons), and whole-transcriptome expression in the same matrix. Caveats:

- Guide capture is sparse: 3.4% of cells in Mouse 2 and <0.5% in Mice 1 and 3; 57% of guide+ cells have a single guide UMI, so the authors' count > 0 rule cannot separate true infection from spill. Mouse 1 and 3 give only tens of cells per guide and will mostly serve as replicates for aggregate statistics, not per-guide models.
- 25% of guide+ cells in Mouse 2 carry more than one guide (authors label them "Multiple" and discard them). In vivo AAV co-infection is plausible, so the multiplet rate is not a clean spillover readout by itself.
- No cell-type annotation or clustering is shipped (`cellTypeID`, `clusterID` all 0); cell types must be derived from expression.
- Mouse 1 (SAW v7, geftool 1.1.2) uses a different segmentation pipeline than Mice 2 and 3 (Deepcell); treat chips as batches.
- Guide UMIs are 0.0044% of total UMI in Mouse 2, so guide-specific spill estimates will have wide intervals; spillover of endogenous marker genes (e.g. Gfap, Olig2, Trem2 into non-expressing neighbours) can be estimated with far more counts on the same cells and used to calibrate.
