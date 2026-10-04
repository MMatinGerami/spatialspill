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

## 2026-10-05 Clonality analysis (exploratory) and its consequence for the recipient policy

reports/audit/perturb_fish/clonality.md, reports/audit/perturb_map/clonality.md.

Perturb-FISH tumour (187,215 cells, 9,374 single-guide calls):
- Same-guide pairs are 4.4 to 5.1 times more frequent than the stratified null in every
  distance bin from 0 to 100 um, with no decay. Same-target-other-guide pairs are only 1.2 to
  1.5 times enriched, different-target pairs 0.8 to 1.1. A flat, guide-specific excess is the
  signature of clonal expansion of transduced tumour cells (the two guides of a target never
  co-occur in a cell, so clones are guide-specific). First-ring-only excess, the misassignment
  signature, is absent.
- 16.8% of assigned cells sit in same-guide connected components of size 2 or more on the
  pruned Delaunay graph (largest 13).
- Unassigned cells adjacent to g cells are shifted towards g's autonomous profile: median
  Pearson correlation 0.72 between (g minus NTC) and (unassigned neighbours of g minus
  unassigned neighbours of NTC) over 30 targets, all above 0.3. Given a 5% guide-call rate and
  clonal growth, most of this is undetected siblings in the "unassigned" pool, not spillover.

Perturb-map: 98% of labelled spots sit in same-phenotype components, as expected for lesion
labels; the unit is the lesion.

Decision (ADR-003): recipients and controls for spillover estimands must carry a confirmed
guide different from g (NTC by default, optionally any other targeting guide). Cells with no
guide call are never recipients in the real-data analyses, because a missing call is
uninformative about perturbation status under clonal growth. This costs power (recipients are
restricted to the few percent of cells with calls) but removes the dominant confound. The
simulator's "clonal" and "misassignment" scenarios test exactly this.

## 2026-10-05 Pre-registration: first real-data spillover analysis (E1)

Datasets: Perturb-FISH tumour (all cells), Perturb-Multi sections 4, 5, 8, 9, 10.
Unit: cell. Strata: sample x cell type. Outcome: log1p size-factor-normalised counts, all
panel genes. Exposure: ring counts of g-cells in bins [0, 15], (15, 30], (30, 60] um
(Perturb-Multi: hepatocyte diameter about 25 um; Perturb-FISH similar), D_max = 60 um.
Recipients: NTC cells (confirmed non-g guide). Controls: NTC cells with no g cell within 60 um.
Clean-control restriction (no perturbed neighbour of any gene) is applied where it leaves at
least 5 controls per stratum, otherwise reported as "not identified".

H1 (calibration, mandatory): treating each NTC guide as a pseudo-target, the fraction of
(guide, ring, cell type, gene) spillover tests with BH q < 0.10 is at most 0.10, and the
fraction with permutation p < 0.05 is within [0.03, 0.07]. Failure means E1 is miscalibrated
on that dataset and no spillover claim is made from E1 there.
H2 (autonomous effects exist and are recoverable): at least 20% of targets in each dataset have
at least one autonomous (gene) hit at q < 0.10 in the pooled cell-type group. This is a
positive control for power, not a scientific hypothesis.
H3 (spillover exists beyond artifacts): after H1 passes, the number of spillover hits at
q < 0.10 for real targets exceeds the number for NTC pseudo-targets scaled by the ratio of
tests, by at least a factor of 3. If not, the dataset's E1 result is "no detectable spillover".
H4 (bleed-through signature): among first-ring spillover hits, the fraction whose gene is also
an autonomous hit for the same target with the same sign is compared with the fraction in
rings 2 and 3; the first-ring excess is the bleed-through estimate. Pre-specified: if the
first-ring fraction exceeds the outer-ring fraction by more than 0.2, the first ring is
reported as artifact-dominated for that dataset.

Analysis script: scripts/run_pipeline.py with configs/e1_perturb_fish.yaml and
configs/e1_perturb_multi.yaml, committed before the results are produced.

## 2026-10-05 Perturb-Multi clonality (exploratory)

reports/audit/perturb_multi/clonality.md (sections 4, 5, 8, 9, 10; 750,000 cells).
Same-guide pair enrichment over the stratified null: 260x within 15 um, 250x at 15 to 30 um,
81x at 30 to 50 um, 23x at 50 to 75 um, 6.8x at 75 to 100 um. Same-target-other-guide pairs
are not enriched (0.7 to 1.3x), so the excess is guide-specific. 47% of assigned cells sit in
same-guide connected components. Unassigned neighbours of g cells carry g's autonomous
signature (median r 0.71 over 60 targets, 95% above 0.3).

Interpretation: a steep, guide-specific decay over about 50 um is consistent with compact
hepatocyte clones (daughter cells stay adjacent) and with local vector spread; both put
undetected g cells next to detected ones. ADR-003 applies. For Perturb-Multi the pool of
confirmed non-g recipients is the 3,828 control cells plus cells with other guides (option
`control_policy: any_other_guide`, to be added to the group definitions).

## 2026-10-05 First pre-registered E1 run on Perturb-FISH, and amendment A1

results/5d30bb4682 (config e1_perturb_fish.yaml, 200 permutations, 98,000 identified tests).
- H1 passes: NTC pseudo-targets give p < 0.05 in 5.2% of tests (autonomous 8.1%, ring 0 7.0%,
  ring 1 3.4%, ring 2 2.4%) and 0% at q < 0.10.
- H2 fails trivially: zero autonomous hits at q < 0.10, and zero spillover hits. The cause is
  resolution, not biology: with 200 permutations the smallest p is 1/201 = 0.005, and BH over
  86,500 tests requires p below 0.1 x k / 86,500 for the k-th smallest p, which 0.005 can reach
  only if about 4,300 tests are all at the floor. The pre-registered design could not have
  detected anything. (Raw rates: 9.2% of autonomous target tests at p < 0.05 versus 8.1% for
  NTC, so the autonomous signal is weak in any case; median 21 treated cells per target.)
- Positivity: ring 0 (0 to 15 um) has a median of 5.5 recipients per target; rings 1 and 2
  have 10.5 and 32.

Amendment A1 (declared before re-running): p-values are permutation-calibrated z-scores. For
each test the permutation null's mean and standard deviation are accumulated, z = (observed
minus null mean) / null sd, and p is two-sided normal. The exact permutation p-value is kept in
the table as `pvalue_perm`. Calibration is re-checked on NTC pseudo-targets (fraction at
p < 0.05 and at q < 0.10) before any hit is interpreted; H1 applies to the z p-values.
H2, H3 and H4 are unchanged.

## 2026-10-05 Amendment A2 (declared before re-running): studentized permutation statistic

The A1 run (results/5d30bb4682, z on raw differences) produced spillover "hits" with
estimate 0.32 and analytic SE 0.33 but p = 1e-72: for sparse exposures (5 recipients in one
stratum) few permutations yield an identified statistic, and the permuted exposed groups have
different sizes from the observed one, so the null sd of the raw difference is not the
sampling sd of the observed statistic (ratio null sd / analytic SE ranged from 0.04 to 70).
Those hits are discarded.

Fix, now the default: the permutation statistic is studentized, t = estimate / analytic SE,
and z = (t_obs - mean t_null) / sd t_null. A test is identified only if at least
max(20, n_perm / 2) permutations produced a valid t. The CI uses the analytic SE rescaled by
sd t_null. On a pure-null toy (4,000 cells, 30 genes, 3 samples, 100 permutations) the
fraction of tests with p < 0.05 is 0.042 (autonomous), 0.047 (ring 0), 0.084 (ring 1), the
exact permutation p gives 0.029, 0.044, 0.078, and 95% CIs cover zero in 93.5% of tests.
H1 to H4 unchanged; H1 is evaluated on these p-values.
