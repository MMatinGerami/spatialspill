# Lab notebook

Dated entries, newest at the bottom. Attempts, failures, pivots and open questions are
recorded here as they happen. Pre-registrations appear under a dated "Pre-registration"
heading before the corresponding real-data analysis is run; anything not pre-registered is
labelled exploratory.

## 2026-10-04 Project start

- Environment: Python 3.11 via uv; squidpy 1.8.2, scanpy 1.11.5, anndata 0.12.19,
  mudata 0.3.10, pertpy 1.0.3, torch 2.14.1 (MPS available), torch_geometric 2.8.0.
- Failure: pertpy 1.0.3 fails to import with statsmodels 0.15.0 (`multipletests` moved).
  Fix: pin statsmodels <0.15 (ADR-002).
- Dataset access verification launched in parallel for Perturb-FISH, Perturb-DBiT
  (GSE319277), Perturb-map, Perturb-Multi, other spatial screens, and supporting resources
  (Replogle via scPerturb, atlases, ligand-receptor databases, embeddings). Results recorded
  in data/README.md when they arrive.
- Compute budget: one Apple Silicon machine, 24 GB unified memory, 15 cores, no CUDA. Heavy
  steps will be subsampled and the downscaling documented.

## 2026-10-04 Dataset access results (Phase 1)

Verified by probing (details, URLs and sizes in data/README.md):

- Perturb-Multi (Cell 2025): best single-cell-resolution screen. 2,206,191 liver cells x 209
  genes, 11 sections, 79,179 guide-assigned singlets over 203 targets including a `control`
  class of 4,161 cells spread over many control guides. Finding: the deposited `raw.X` is
  log1p(normalised), not counts. The per-cell size factor is recoverable from the smallest
  non-zero value (count 1), and dividing gives integers to 1e-13 with target_sum 93. Counts are
  recovered exactly in the loader (unit test added). Downscaling: sections 4, 5, 8, 9, 10 (the
  five with more than 10,000 assigned cells) and at most 150,000 cells per section for the
  audit; the estimator runs will state their own subsets.
- Perturb-DBiT (Nat Biotechnol 2026): pixel-level (20 to 50 um microfluidic pixels), 6 of 14
  samples have both expression and sgRNA tables in GEO. 23% of pixels carry a dominant guide,
  so the "clean control" restriction (no perturbed neighbour within D_max) leaves no
  controls; the DBiT analysis uses all unperturbed pixels as controls and treats other-guide
  exposure as a covariate (E2) instead. Only 4 NTC pixels in mTSG.LV.2: NTC calibration is
  not possible there; Ep.Tu.1 has NTC guides but no expression matrix in GEO.
- Perturb-map (Cell 2022): Visium spot level with lesion Pro-Code phenotypes for a subset of
  knockouts; per-cell imaging tables on BioImage Archive S-BIAD267. Clonal lesion design;
  lesion is the unit.
- Perturb-FISH (Cell 2025): tables on Brain Image Library. The in vivo tumour xenograft has
  coordinates, counts and a 77-column perturbation design linked by row order; THP-1 has no
  per-cell coordinates deposited (only neighbour counts), so the in vitro intercellular
  analysis of the paper cannot be reproduced at the cell level without re-segmentation.
- Spatial Perturb-seq (Nat Commun 2026): Stereo-seq cellbin GEF downloaded (2.1 GB); guide
  features inside the GEF still to be verified.
- Not public: PerturbView, PerturbSpace, SPACE.
- Infrastructure failures: GEO's HTTPS front end returned 403 intermittently; switched to
  ftp://. The first CI run failed on a ruff rule (RUF007) fixed in the next commit; the
  second needed the smoke pipeline that did not exist yet.

Open question: at 3.6% assigned cells, how many Perturb-Multi recipients have exactly one
perturbed neighbour of a given gene within 30 um? The audit and the first E1 run answer this
(positivity table).

## 2026-10-04 Perturb-Multi audit: same-guide cells are strongly clustered

Audit of sections 4, 5, 8, 9, 10 (750,000 cells after capping at 150,000 per section;
reports/audit/perturb_multi/). Density is balanced between perturbed, NTC and unassigned cells
(5.53, 5.54, 5.48 neighbours within 30 um), so delivery is not density-biased. But the number
of graph edges joining two cells with the same target is 24,856 on the pruned Delaunay graph
against a stratified-permutation null of 181 (sd 15): z about 1,700, and similar on the radius
and kNN graphs. Same-target cells are about 140 times more often adjacent than chance.

Candidate explanations, to be separated before any spillover claim:
1. Clonal expansion of transduced hepatocytes (biology; would make neighbours siblings).
2. Local spread of the viral vector along sinusoids (delivery; neighbours independently
   perturbed with the same guide; still breaks A3 within strata).
3. Barcode misassignment: RCA amplicons from one cell detected in a neighbour (measurement;
   the "neighbour" is not perturbed at all, and its "autonomous" effect is bleed-through).

Consequence: a recipient's exposure to gene g is strongly correlated with being itself a
(possibly unlabelled) g-cell. Spillover estimates that do not model this will be inflated by
the autonomous effect. Planned analysis (exploratory, then pre-registered tests):
distance-decay of same-guide pairs vs same-gene-different-guide pairs (clonality and delivery
predict same-guide; misassignment predicts same-guide too but confined to the first ring), and
whether unassigned cells adjacent to g-cells show g's autonomous expression signature.
