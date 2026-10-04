# Perturb-FISH tumour xenograft: raw table inspection

Dataset: Binan et al., Cell 2025, Brain Image Library doi:10.35077/ace-gem-get, folder
`extras/tumors/processed/finaltables/` (one FFPE section of a humanised NSG mouse tumour
xenograft, MERSCOPE 500-gene immuno-oncology panel, 35-target TLR/NF-kB CRISPR library).
Author code: github.com/lbinan/Perturb-FISH (`tumorpreprocessing/`, `tumoranalysis/`).
All numbers below were produced by the scripts named in each section (kept in the session
scratchpad, `forensic1.py` to `forensic5.py` and `run_loader.py`), run on 2026-10-04 with
`uv run python`.

## Files

| File | Shape | Header | Content |
|---|---|---|---|
| `merfishcounttable.csv` | 187,215 x 553 | no | col 1 cell index 1..N; col 2 area (mosaic px); col 3 MATLAB linear index of centroid; cols 4..503 gene counts (500 genes); cols 504..553 blank barcodes (50) |
| `coordinates.csv` | 187,215 x 2 | no | centroid x (mosaic column), y (mosaic row), mosaic pixels |
| `allcellsPerturbationTable.csv` | 187,215 x 77 | no | binary guide detections per cell (values 0/1 only) |
| `tumorMerfish.csv` | 9,432 x 500 | yes (gene names) | counts of the authors' perturbed, QC-passing cells |
| `tumorpooledperturbations.csv` | 36 x 9,432 | yes (`Cell_1..`) | pooled design (35 targets + `Control`) for those cells, values 0/1/2 |
| `withimmuneneighborMerfish.csv` / `withoutimmuneneighborMerfish.csv` | 3,593 x 500 / 5,839 x 500 | yes | the 9,432 perturbed cells split by having a T-cell neighbour |
| `withneighborperturbation.csv` / `withoutnotimmuneneighborperturbations.csv` | 36 x 3,593 / 36 x 5,839 | yes | matching designs |
| `perturbed_tcells_merfish.csv` | 14,232 x 113 | yes | counts (113 genes with mean > 0.3 in T cells) for the authors' T cells |
| `perturbed_tcells_perturbation.csv` | 19 x 14,232 | yes | for each T cell, perturbation of neighbouring tumour cells (18 targets + `Control`; 12,966 cells are control-only, i.e. no perturbed neighbour) |
| `tumorcells_LFCs.csv`, `*_qvals.csv`, `split*`, `tcells_*` | 500 x 36, 113 x 19, ... | yes | FR-Perturb effect tables; column headers give the 35 target names + `Control` |
| `tumor_xenograft/` | | | byte-identical copies of the above (same sizes) |
| `../sample1/MERFISH/DeriveEntityMetadataTask/region_0/cell_metadata.csv` | 419,567 x 16 | yes | Vizgen segmentation: fov, volume, center_x/y (um), ... (fetched to the scratchpad; the project download script fetches it last) |

The MERSCOPE `experiment.json` reports 1,020 FOVs and 419,567 Vizgen cells (74,482 with no
transcripts). The authors re-segmented the DAPI mosaic themselves (watershed, `tumormerfishcounttable.m`),
which is why their table has 187,215 cells.

## How the tables link

### coordinates.csv is row-aligned with merfishcounttable.csv

`tumormerfishcounttable.m` stores `sub2ind(size(mymask), Centroid(2), Centroid(1))` in column 3,
and `preparetablesforPython.m` notes `size(mymask)=169973 102389`. With row = y and
column = x, `idx = y + (x - 1) * 169973` holds for every row (forensic1.py):

```
coordinates shape (187215, 2)
linkage residual |idx - (y+(x-1)*169973)|: max 5.9e-05 (all rows)
col1 == 1..N: True
```

The alternative decompositions (x as row, or 102,389 as the leading dimension) give residuals
of 1e8 to 1e10. So coordinates.csv row i is the centroid of count-table row i, in mosaic pixels
(x range 37 to 91,986; y range 92,171 to 153,204).

### Gene names: tumorMerfish.csv header = count-table columns 4..503

`preparetablesforPython.m` uses `counts = merfishcounttable(:, 4:end); counts = counts(:, 1:500)`.
Columns 504..553 are the 50 blank barcodes: their per-column sums are low (median 3,316, max
11,605, total 194,144) against a median of 12,018 per gene column (min 1,074), total 38,019,784.
The header of `tumorMerfish.csv` has exactly 500 names (PDK4, CCL26, ..., NOTCH1), identical to
the row names of `tumorcells_LFCs.csv`. Every one of the 9,432 rows of `tumorMerfish.csv` matches
exactly one row of `merfishcounttable[:, 4:503]` (forensic3.py):

```
tumorMerfish rows 9432 matched 9432 ambiguous 0 missing 0
matched indices monotone increasing? True
matched cells: tot range 41 798 area range 2105.0 29995.0
```

The matched cells lie inside the authors' filter `totalcounts > 40 & < 800 & area > 2000 & < 30000`,
which keeps 153,705 of 187,215 cells, the exact `n_obs` printed in `tumoranalysis/clusterCells.ipynb`.
Since the match used all 500 columns in header order, the header order is the count-table order.

### allcellsPerturbationTable.csv: 77 guide columns

Column sums (1-based):

```
1 13 | 2 21 | 3 40 | 4 15 | 5 64 | 6 14 | 7 1046 | 8 132 | 9 22 | 10 301 | 11 108 | 12 320 |
13 189 | 14 126 | 15 256 | 16 39 | 17 19 | 18 20 | 19 54 | 20 61 | 21 13 | 22 25 | 23 7 | 24 54 |
25 275 | 26 5 | 27 258 | 28 11 | 29 2 | 30 991 | 31 117 | 32 208 | 33 7 | 34 70 | 35 25 | 36 38 |
37 22 | 38 6 | 39 27 | 40 104 | 41 20 | 42 102 | 43 182 | 44 17 | 45 40 | 46 176 | 47 60 | 48 1 |
49 115 | 50 14 | 51 123 | 52 79 | 53 21 | 54 14 | 55 26 | 56 109 | 57 55 | 58 1 | 59 1 | 60 24 |
61 74 | 62 43 | 63 1 | 64 58 | 65 16 | 66 10 | 67 55 | 68 1147 | 69 162 | 70 42 | 71 12 | 72 2822 |
73 1039 | 74 1398 | 75 0 | 76 12 | 77 1
```

`preparetablesforPython.m` gives the layout:

```matlab
corrected(:,7)=0;corrected(:,30)=0;corrected(:,68)=0;  % noisy barcodes
corrected(:,71)=corrected(:,71)+corrected(:,72)+corrected(:,73)+corrected(:,74);
corrected=corrected(:,1:71);
for i=1:35
    pooled(:,i)=corrected(:,2*i-1)+corrected(:,2*i);   % 2 guides per target
end
pooled(:,36)=corrected(:,71);                           % control
```

So columns 2i-1 and 2i are the two guides of target i, columns 71..74 are control guides,
columns 75..77 (sums 0, 12, 1) are unused. Target order i = 1..35 is given by the index map in
`clusterCells.ipynb` and by the row names of `tumorpooledperturbations.csv`:

CD14, CHUK, IKBKB, IRAK1, IRAK4, IRF3, IRF5, IRF7, JUN, LBP, LY96, MAP2K1, MAP2K2, MAP2K3, MAP2K4,
MAP2K6, MAP2K7, MAP3K7, MAPK14, MYD88, NFKB1, NFKBIA, PELI1, PIK3CA, RELA, RIPK1, TAB1, TAB2, TBK1,
TICAM1, TIRAP, TLR4, TRADD, TRAF6, TRAM1, Control.

Independent check of the mapping: for the 9,432 published perturbed cells, the raw column set in
`allcellsPerturbationTable.csv` was compared with the pooled design row (forensic3.py). Columns map
to the expected target in 80 to 100 percent of cells (e.g. col 8 -> IRAK1 0.97, col 10 -> IRAK4
1.00, col 25 -> MAP2K2 0.99, col 46 -> PELI1 1.00, cols 72/73/74 -> Control 1.00), while the three
"noisy" columns do not (col 7 -> Control 0.22, col 30 -> Control 0.39, col 68 -> TAB2 0.26),
confirming that they are not IRAK4_g1, MAP2K4_g2, TRAF6_g2 detections. Over all cells with a
single raw guide and a single published target, agreement is 98.4 percent (6,508 cells).
A side observation: the two guides of the same target never co-occur in a cell (0 for all 35
pairs, against a mean of 0.16 co-occurrences for arbitrary column pairs).

### The published design is not identical to the raw table

The GitHub rule applied to `allcellsPerturbationTable.csv` yields 8,603 perturbed QC-passing cells,
but `tumorMerfish.csv` has 9,432, and 2,547 of those have no guide at all in the raw table
(forensic2.py tried every combination of zeroed/control columns and two count thresholds; none
reproduces 9,432). The published design therefore comes from a later guide decoding ("more cells")
than the deposited `allcellsPerturbationTable.csv`. The loader uses the raw table for `guide` and
`target` (guide-level, all cells) and stores the authors' pooled call as `obs["target_published"]`
for the 9,432 matched cells.

### Author labels recovered by exact expression matching (forensic4b.py)

```
perturbed_tcells_merfish rows 14232 unique match 14232 ambiguous 0 missing 0
with immune nb rows 3593 matched 3593 amb 0 miss 0
without immune nb rows 5839 matched 5839 amb 0 miss 0
union with/without == tumorMerfish matched set: True overlap 0
T cells overlapping perturbed tumour set: 264
```

The 14,232 T cells are the authors' leiden cluster 5 (TRAC, CD2, PTPRC, GZMB, CD8A markers in
`clusterCells.ipynb`), sub-clustered and with "bad cells" removed in MATLAB. They are labelled
`cell_type = "T cell"` in the loader; all other cells are `"unknown"` (no clustering is done here).

### Units: mosaic pixels to micrometres (forensic4a.py, forensic5.py)

No `micron_to_mosaic_pixel_transform.csv` is deposited for the author mosaic. Comparing the
0.1 to 99.9 percentile spans of the author coordinates with those of Vizgen `cell_metadata.csv`
centres (same section):

```
px span x,y 89890 60311   um span x,y 9791 6369
ratios px/um: x/x 9.18  y/y 9.47   (x/y 14.1, y/x 6.2: axes are not swapped)
```

Both ratios bracket the MERSCOPE mosaic pixel of 0.108 um (9.26 px/um), so the loader uses
`um_per_px = 0.108` (2 to 3 percent uncertainty). With this scale the author median area of
11,941 px^2 becomes 139 um^2 (Vizgen median volume 98.8 um^3, in a different segmentation).
A per-cell match to Vizgen cells was attempted (2D histogram cross-correlation over rotations
-3 to +3 degrees, then nearest neighbour): median nearest-neighbour distance 8.2 um, 25 percent
of cells within 5 um. That is too loose to transfer `fov` or `volume`, so `batch = sample` and
`area` is the author pixel area times 0.108^2.

## Loader output (run_loader.py)

```
load time s 3.8
AnnData object with n_obs x n_vars = 187215 x 500
validate -> []
spatial range um x 4.0 .. 9934.5   y 9954.5 .. 16546.1
```

Guides per cell (raw table after zeroing columns 7, 30, 68, 75, 76, 77):

| n_guides | cells | percent |
|---|---|---|
| 0 | 177,533 | 94.83 |
| 1 | 9,374 | 5.01 |
| 2 | 298 | 0.16 |
| 3 | 10 | 0.01 |

Cells with 2 or more guides get `guide = "multi"`, `target = "none"`, `is_perturbed = False`.
`is_ntc` 4,987 cells; `is_perturbed` 4,387 cells; 35 distinct targets.

Counts per target (single-guide cells):

```
NTC 4987 | IRF3 402 | MAP2K6 309 | IRAK4 303 | IRF7 280 | IRF5 277 | MAP2K2 266 | MAP2K3 239
PELI1 210 | RIPK1 186 | NFKBIA 184 | TRAM1 182 | IRAK1 127 | MYD88 125 | TAB2 124 | RELA 118
TIRAP 112 | NFKB1 109 | LBP 104 | IKBKB 72 | MAP2K7 71 | PIK3CA 61 | MAP2K1 59 | MAP3K7 57
TLR4 55 | TRAF6 53 | TBK1 51 | CHUK 51 | JUN 37 | LY96 34 | TAB1 32 | CD14 27 | MAPK14 26
TRADD 22 | TICAM1 21 | MAP2K4 1
```

Other obs: `author_qc_pass` 153,705; `cell_type` T cell 14,232 / unknown 172,983;
`target_published` set for 9,432 cells (NTC 4,264, multi 261); `has_immune_neighbor`
yes 3,593 / no 5,839 / unknown 177,783; area median 139 um^2 (5th to 95th percentile 53 to 318);
total_counts median 160; blank_counts median 1.

## Uncertain points

- Pixel size 0.108 um/px is inferred from span ratios (9.18 to 9.47 px/um), not from a deposited
  transform; absolute positions are in the author mosaic frame, not Vizgen's.
- `allcellsPerturbationTable.csv` is an earlier decoding than the one behind the published
  9,432-cell design (2,547 published cells have no raw guide). Results from `target` and
  `target_published` will differ slightly.
- Guide-level names (`CD14_g1`, ...) are positional (2i-1, 2i); the deposited data give no sgRNA
  sequences, so g1/g2 cannot be tied to specific protospacers.
- MAP2K4 has 1 usable cell (its second guide, column 30, was zeroed as noisy); the published
  design also has only 2 MAP2K4 cells.
- 264 of the authors' T cells carry a guide call (110 targeting, 77 control in the raw table).
  The authors treat guide-carrying cells as tumour cells; these may be doublets or segmentation
  spill. They are labelled "T cell" here because that is the explicit author label.
- `area` is a 2D segmentation area from the author watershed mask, not Vizgen's 3D volume.
