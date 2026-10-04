# Data inventory

Verified 2026-10-04 by probing each resource (HTTP HEAD, GEO FTP listings, Crossref, repository
APIs). "Verified" means the item appeared in a tool result. md5 sums and cell/perturbation counts
are filled in by `scripts/download_data.py` after download (written to `data/raw/MANIFEST.md`).
Nothing in `data/raw` is tracked by git.

## Core spatial perturbation screens

### Perturb-Multi (Zhuang, Weissman and Allen labs), Cell 2025. Accessible: yes
- Paper: Saunders RA et al., Perturb-Multimodal: a platform for pooled genetic screens with
  imaging and sequencing in intact mammalian tissue. Cell 2025. doi:10.1016/j.cell.2025.05.022
- Data: Hugging Face dataset `xingjiepan/PerturbMulti` (public, not gated, CC BY 4.0).
  - `RNA_scaled_crispr_screen_20240615.h5ad` (14.2 GB): 2,206,191 cells x 209 MERFISH genes,
    `raw` layer present; obs includes `x`, `y`, `global_x`, `global_y`, `fov`, `batch` (11),
    `area`, `cell_type` (9), `singlet_gene` (204 target genes), `singlet_name` (456 guides).
  - `protein_intensities_crispr_screen_20240615.h5ad` (31 MB): 99,294 guide-assigned cells x 18
    protein channels with `x`, `y`, `perturbation`.
  - Per-cell images in 42 tar chunks (~375 GB): not used.
- Sequencing arm: GEO GSE275483 (10x Perturb-seq of liver), not spatial.
- Code: github.com/weallen/InVivoMultimodalPerturbation (no licence file).
- Role here: primary single-cell-resolution development dataset (mouse liver, in vivo).
  Downscaling: the full h5ad is read with h5py; only cells in fields of view that contain at
  least one guide-assigned cell are kept.

### Perturb-DBiT (Fan lab), Nat Biotechnol 2026. Accessible: yes
- Paper: Baysoy A et al., Large-scale, spatially resolved panoramic CRISPR screening in native
  tissue environments using Perturb-DBiT. Nat Biotechnol 2026. doi:10.1038/s41587-026-03127-y
  (CC BY 4.0). Preprint doi:10.1101/2024.11.18.624106. The Nat Biotechnol 2026 piece "Spatial
  CRISPR screens map total RNA in tissue" (doi:10.1038/s41587-026-03126-z) is a News and Views
  on this paper, not a separate dataset.
- GEO GSE319277 (public 2026-03-17; 14 samples; RAW tar 78 MB). Per sample: `position_*.txt.gz`
  (pixel grid indices AxB, 50x50 or 100x100 microfluidic channels of 20 or 50 um),
  `sgCounts_*.txt.gz` (columns BA, BB, grna, count), ROI PNG; dense pixel x gene TSV for 8 of 14
  samples (Bru.LC.4, Bru.LC.5, Bru.LC.6_allRNA, mTSG.LV.2, mTSG.LV.3, Bru.LC.8, Brie.LU.1,
  Bru.LC.1.33.3). No expression matrix in GEO for NTC.SP.1, Ep.Tu.1, Bru.LC.1/2/3/7 (FASTQ only,
  SRA PRJNA1422738).
- GEO GSE319123: one Xenium Prime 5K section (HT29 lung metastasis, FFPE) with a custom 27-sgRNA
  probe panel. Files `cells.parquet.gz` (6 MB), `cell_feature_matrix.h5` (14 MB),
  `cell_boundaries.parquet.gz` (30 MB), `transcripts.parquet.gz` (180 MB), `morphology.ome.tif.gz`
  (13 GB, not used).
- Code: github.com/abaysoy/Perturb_DBiT (no licence file).
- Role here: sequencing-based technology. Unit of analysis is the pixel (not the cell); the
  Xenium section gives a single-cell-resolution companion with 27 guides.

### Perturb-map (Merad and Brown labs), Cell 2022. Accessible: yes
- Paper: Dhainaut M et al., Spatial CRISPR genomics identifies regulators of the tumor
  microenvironment. Cell 2022. doi:10.1016/j.cell.2022.02.015. PMID 35290801.
- GEO GSE193460 (RAW tar 150 MB): 4 Visium sections of KP lung tumours (3 mice), Space Ranger
  outputs plus `spot_annotation.csv.gz` with `kmeans` (normal/tumour), `leiden_clusters` and
  `phenotypes` (Pro-Code identity for highlighted lesions: Tgfbr2_1, Jak2_1, KP clusters).
- BioImage Archive S-BIAD267 (15.7 GB, 164 files): per-cell QuPath intensity tables
  (`*_segData.csv`, centroid x/y, nucleus area, marker and Pro-Code epitope intensities;
  Lymphoid and Myeloid panels) and curated tumour boundary JSON for 41 lung lobes. Pro-Code
  debarcoding is done by `perturb_map_processing_functions.R` in the code repo.
- Code: github.com/srose89/Perturb-map (MIT). No scRNA-seq deposit exists for this paper.
- Role here: clonal-lesion design; unit of analysis is the lesion; validation ground for C6.

### Spatial Perturb-seq (Chew lab), Nat Commun 2026. Accessible: yes (with caveat)
- Paper: Shen K et al., Spatial perturb-seq: single-cell functional genomics within intact tissue
  architecture. Nat Commun 2026. doi:10.1038/s41467-026-69677-6. Preprint doi:10.1101/2024.12.19.628843.
- GEO GSE274447: 3 Stereo-seq cellbin GEF files (mouse hippocampus, AAV library of 18 targets),
  RAW tar 2.1 GB. Whether guide barcodes are stored as features inside the GEF is unverified
  until download. Code: github.com/kimberle9/spatialperturbseq, Zenodo doi:10.5281/zenodo.17959756.
- Role here: candidate fourth technology if the GEF carries guide features.

### Combinatorial functional genomics + spatial transcriptomics (Gerstung and Tschaharganeh labs), Nat Biomed Eng 2025. Accessible: yes
- Paper: Breinig M et al. doi:10.1038/s41551-025-01437-1 (CC BY 4.0). Zenodo record 10986436
  (`data_tidy.zip`, 993 MB, CC BY 4.0). Visium spot-level, genotype per spot. Not used in the
  first pass (spot-level unit, combinatorial genotypes); listed for completeness.

### Perturb-FISH (Binan et al.), Cell 2025
- Paper: doi:10.1016/j.cell.2025.02.012. Access verification in progress (see below).

### Not public (as of 2026-10-04)
- PerturbView (doi:10.1038/s41587-024-02391-0): images "will be available on Image Data Resource";
  no IDR record found.
- PerturbSpace (doi:10.64898/2026.05.25.727765): data upon publication.
- SPACE (doi:10.1101/2025.09.14.675819): no accession found.

## Supporting data
(Filled in when verification completes.)
