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
  RAW tar 2.1 GB. Guide barcodes are stored as 18 `sgrna_<target>` features in the gene table
  of each GEF (verified after download, reports/audit/spatial_perturbseq_inspection.md). Code: github.com/kimberle9/spatialperturbseq, Zenodo doi:10.5281/zenodo.17959756.
- Role here: fourth technology (Stereo-seq cell bins); exploratory E2 and clonality analyses.

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

### Perturb-FISH (Binan et al.), Cell 2025. Accessible: yes (Brain Image Library)
- Paper: doi:10.1016/j.cell.2025.02.012 (CC BY). PMID 40081369. Preprint doi:10.1101/2023.11.30.569494.
- Data: Brain Image Library doi:10.35077/ace-gem-get (CC BY 4.0), anonymous HTTPS directory
  `https://download.brainimagelibrary.org/0c/bd/0cbd479c521afff9/`. About 26 TB in total (raw
  images); the processed "finaltables" are small. GEO GSE221321 is the matched compressed
  Perturb-seq (scRNA-seq), not spatial. No h5ad is provided; tables are linked by row order.
  - Tumour xenograft (human melanoma, NF-kB pathway knockouts, humanised NSG mice, FFPE, 500-gene
    immuno-oncology MERFISH panel): `coordinates.csv` (187,214 cells, x/y in microns),
    `merfishcounttable.csv` (212 MB), `allcellsPerturbationTable.csv` (77 binary columns), design
    matrices for pooled/neighbour analyses, FR-Perturb effect tables, and MERSCOPE
    `cell_metadata.csv` (113 MB; volume, centre, fov). This is the in vivo dataset with
    published intercellular effects and anchors C3/C4.
  - THP-1 (LPS-stimulated, 74 gRNAs / 35 targets, 129 genes): count tables for 7,088 and
    8,097 cells with target-level design matrices; per-cell coordinates are not deposited
    (only `numberofneighbors` for sample 2); would need re-segmentation from raw images.
  - Astrocytes (iPSC-derived CRISPRi, 127 genes, 277-gene panel): 14,926 cells, design matrix,
    MERSCOPE cell_metadata with centroids; calcium clusters.
- Code: github.com/lbinan/Perturb-FISH (GPL-3.0), github.com/douglasyao/FR-Perturb.
- Download: `scripts/download_perturb_fish.sh` (final tables + cell metadata only).
- Caveat: guide-level calls are not in the tables (target level only); linkage of counts,
  coordinates and design is by row order as documented in the README and the MATLAB scripts.

## Supporting data (verified 2026-10-04)

| Resource | Location | Size | Licence | Use |
|---|---|---|---|---|
| Replogle 2022 pseudobulk (K562 essential, K562 gwps, RPE1; raw and normalised) | figshare doi:10.25452/figshare.plus.20029387, files 35773070 / 35780870 / 35774443 / 35773217 / 35775581 / 35775512 | 80 to 375 MB each | CC BY 4.0 | priors on autonomous effects; perturbation-derived gene embeddings |
| Replogle 2022 single cell via scPerturb | Zenodo doi:10.5281/zenodo.13350497 (`ReplogleWeissman2022_*.h5ad`) | 1.2 to 8.8 GB | CC BY 4.0 | not needed if pseudobulk suffices |
| Allen Brain Cell Atlas MERFISH (Yao 2023, 500 genes), per-section h5ad | `s3://allen-brain-cell-atlas/expression_matrices/MERFISH-C57BL6J-638850-sections/20240330/` | 23 to 261 MB per section | Allen terms (CC BY 4.0) | empirical null spatial correlation (brain) |
| 10x Xenium public: mouse brain (v1, 5K), human lung cancer FFPE (v1, 5K) | cf.10xgenomics.com/samples/xenium/... `cell_feature_matrix.h5`, `cells.parquet` | 3 to 104 MB (matrices) | 10x dataset terms (CC BY 4.0, not re-verified) | empirical null (lung, brain) |
| Vizgen MERFISH mouse liver map | registration wall; GCS 401/403 | n/a | not verified | not usable without registration; Perturb-Multi's own unassigned cells serve as liver null |
| OmniPath intercell + ligrecextra | omnipathdb.org (TSV) | 2 to 7 MB | academic | ligand-receptor enrichment |
| CellChatDB human/mouse | github.com/jinworks/CellChat `data/*.rda` | 1.3 to 1.5 MB | GPL-3 | ligand-receptor enrichment |
| CellPhoneDB v5 | github.com/ventolab/cellphonedb-data (`interaction_input.csv` etc.) | 0.5 MB | not stated | ligand-receptor enrichment |
| NicheNet v2 ligand-target matrices (human, mouse) | Zenodo doi:10.5281/zenodo.7074291 | 191 to 262 MB | CC BY 4.0 | sender-receiver prior |
| ESM2 checkpoints | huggingface.co/facebook/esm2_t12_35M_UR50D (136 MB), esm2_t6_8M (31 MB) | | MIT | protein embeddings for C5 |
| GO annotations and ontology | current.geneontology.org `goa_human.gaf.gz` (15 MB), `go-basic.obo` (32 MB) | | CC BY 4.0 | GO embeddings for C5 |
| DepMap 24Q4 CRISPRGeneEffect | figshare article 27993248, file 51064667 | 429 MB | CC BY 4.0 | target plausibility (newer releases are portal-only behind a bot check) |
| Open Targets Platform 26.09 | GraphQL api.platform.opentargets.org/api/v4/graphql; bulk ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/ | | CC0 | target plausibility |
| scPerturBench metrics | github.com/bm2-lab/scPerturBench `Perturbation_generalization/calPerformance_genetic.py` | | GPL-3 | C5 metric definitions (pertpy Distance on top DEGs) |
