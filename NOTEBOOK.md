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

## 2026-10-05 Results of the A2 run on Perturb-FISH (results/5d30bb4682) and amendment A3

Run: configs/e1_perturb_fish.yaml, studentized E1, 200 permutations, 95,000 identified tests.
- H1 (NTC calibration): at q < 0.10, 1.0% of NTC tests are called (criterion: at most 10%,
  passes). At p < 0.05: 8.1% autonomous, 13.4% ring 0, 6.1% ring 1, 7.3% ring 2 (criterion
  [0.03, 0.07]: fails for ring 0 and marginally for autonomous and ring 2). E1 is
  anticonservative, most strongly where exposure is sparsest (ring 0: median 7.5 recipients).
- H2: 32 of 32 targets have at least one autonomous hit at q < 0.10 (1,052 hits). Passes,
  but with the same anticonservativeness caveat.
- H3: 637 spillover hits against an NTC rate that predicts about 84,000 x 0.009 = 760 false
  calls; the criterion (3x the NTC-scaled count) fails. No spillover claim from E1.
- Diagnosis: the strongest ring-0 "hits" are MAP2K6 on BST2, AKT1, CD40, ICAM1, all about
  -0.5 with exactly 5 recipients. Five NTC cells adjacent to one MAP2K6 clone share one
  niche; their joint deviation reflects local expression autocorrelation, which the
  stratified permutation null (strata = sample x cell type, one sample, one cell type for
  tumour cells) does not model: it compares them with 5 random NTC cells anywhere in the
  section. This is assumption A4 failing through spatial autocorrelation, and it is the
  general weakness of E1 with global strata.

Amendment A3 (declared before running): strata become sample x cell type x spatial tile,
with square tiles of 250 um (about 10 cell diameters, chosen so that a tile holds a clone
and its ring recipients together with local NTC controls). Differences and permutations are
then local. Positivity will drop (tiles with fewer than 5 recipients or controls are
unidentified) and this is reported. The tile size sensitivity (150, 250, 400 um) is run and
reported; the pre-specified primary is 250 um. H1 is evaluated again under A3; H3 and H4 are
evaluated only if H1 passes.

Benchmark note (results/bench_small_a2): on the simulator, E1's 95% intervals cover zero for
94 to 95% of null spillover pairs (ring 2) but only 89 to 91% of null autonomous pairs, with
FDP 0.18 to 0.24 at q < 0.10 for autonomous effects. The autonomous "null" pairs are not null
on the analysis scale: a perturbation that changes 10% of genes shifts the cell's total
counts, and size-factor normalisation moves every other gene (compositional bias). The
benchmark now also scores NTC pseudo-targets, which are exact nulls, and the manuscript will
report autonomous calibration on them. The compositional effect is itself a finding for
panel-based assays (209 to 500 genes), where one gene can be a large share of the total.

## 2026-10-05 E4 (GNN) implemented; its uncertainty is not usable

E4 (src/spatialspill/estimators/e4_gnn.py): 2-layer message passing on the pruned Delaunay
graph with own-guide one-hots as node features; spillover from single-neighbour interventional
counterfactuals. On the simulator it recovers the direction and ranking of planted ring-0
spillover (r = 0.91 over 20 target-gene pairs, 4,000 cells) but shrinks magnitudes (slope
0.09), and its recipient-bootstrap standard errors ignore the variance of the fitted function:
on a pure null 73% of identified tests have p < 0.05 (E2: 4.3%). E4 therefore enters the
benchmark for ranking metrics (AUROC) only; its p-values and intervals are not reported as
inference. Options noted for later: retraining ensembles, or calibrating against the NTC
pseudo-target null.

## 2026-10-05 Tile strata on Perturb-FISH (amendment A3), 250 um

results/8c106a7091: 8,000 identified tests (from 95,000), NTC p < 0.05 fraction 0.036,
q < 0.10 fraction 0; 0 autonomous hits, 4 spillover hits. Local strata remove the
autocorrelation-driven false positives and most of the power with them. 150 and 400 um runs
pending; whatever they show, E1 with global strata is not an acceptable spillover test on a
single-section dataset with clonal structure, and the manuscript will say so.

## 2026-10-05 E4 benchmark and pause

results/results/0603f268a6 (10,000 cells on Perturb-FISH geometry, 15 targets, 30% assigned, strong planted
spillover, any_other_guide recipients, 2 reps): for spillover detection E4 ranks far better
than E1 and E2 (AUROC 0.88 to 0.92 vs 0.53 to 0.63 for E1/E2 at this size) but its p-values
are useless (null FPR above 0.9), while E1's NTC pseudo-targets are calibrated (5.3 to 6.1%
at p < 0.05, 94 to 95% coverage) and E2 is slightly anticonservative in ring 0 (12.9%).
Conclusion so far: ranking and inference need different tools; E4 for ranking, E1 with local
strata for calibrated claims. The 40,000-cell power benchmark with E0 to E3 (results/bench_power2.log),
the Perturb-FISH tile sensitivity runs (150 and 400 um) and the Perturb-Multi E1 runs
(global and 400 um tiles) were still running in tmux when work paused on 2026-10-05 at 01:30.

## 2026-10-05 07:37 Overnight results

Perturb-FISH tile sensitivity (E1, NTC recipients):
- 150 um (results/a16b3fdbc3): 4,000 identified tests, NTC p < 0.05 fraction 0.025, 0 hits.
- 250 um (results/8c106a7091): 8,000 tests, 0.036, 0 autonomous, 4 spillover hits.
- 400 um (results/055d929ee0): 15,000 tests, 0.038, 0 autonomous, 6 spillover hits.
- global (results/5d30bb4682): 95,000 tests, 0.083 (ring 0: 0.134), 1,052 autonomous, 637 spillover.
With local nulls the autonomous hits vanish too: a perturbed clone compared with NTC cells
elsewhere in the section differs by niche as much as by genotype. In a clonal tumour screen
"autonomous effect" estimates that are not locally controlled are confounded, not only the
spillover ones. H2 as pre-registered therefore fails under A3; the honest statement is that
E1 cannot separate genotype from niche in this dataset at this sample size.

Perturb-Multi (sections 4, 5, 8, 9, 10; 5 sections, 9 cell types, 50 NTC guides):
- global strata (results/53c85fb53d): 66,044 tests, NTC p < 0.05 fraction 0.111, q < 0.10
  fraction 0.025; 1,579 autonomous and 94 spillover hits. Anticonservative; H1 fails.
- 400 um tiles (results/7c7029afd1): 0 identified tests. With 3,828 NTC cells over 5 sections
  x 9 cell types x tiles, no stratum reaches 5 recipients and 5 controls. Tile strata are not
  viable where cell types are many; the spatial adjustment has to be a covariate, not a stratum.

Power benchmark, 40,000 cells on Perturb-FISH geometry, any_other_guide recipients, no clean
restriction (results/a0b119699d): for spillover the neighbour t-test (the published approach)
has the highest apparent power (0.6 to 0.7 at ring 0) but FDP 0.27 to 0.59 and 10% NTC false
positives under clonal assignment; E1 and E2 keep NTC false positives at 4 to 7% except E1
in the clonal scenario at ring 1 (14 to 15%), where E2 stays at 6 to 7% (its density and area
covariates absorb part of the clustering); E3 is anticonservative at ring 0 (13 to 16%),
likely from the ridge outcome model and small exposed groups. Spillover power at q < 0.10 is
0.3 to 0.5 in rings 1 and 2 for E1/E2 with strong planted effects (LFC sd 1.5 on 20% of genes),
AUROC 0.8 to 0.87. Autonomous NTC "false positives" of 0.2 to 0.3 under any_other_guide are a
property of that policy: controls then include perturbed cells with their own autonomous
effects, so the NTC-vs-control contrast is not null; under the ntc policy (bench_small_a2) the
NTC rate is nominal. This is recorded as a policy caveat, not an estimator failure.

Decision: implement the pre-registered spatial random effect in E2 as a per-sample radial
basis (k-means centres on coordinates, Gaussian bumps) entering the regression as covariates,
so controls stay available everywhere while niche variation is absorbed. Pre-registration for
its real-data use follows below once the simulator shows it is calibrated under clonal
assignment.

## 2026-10-05 07:42 Pre-registration: E2 with spatial basis on real data

E2 (regression adjustment) with strata sample x cell type, covariates local density, log area,
edge distance, and a per-sample Gaussian radial basis of 40 k-means centres (bandwidth =
median centre spacing), recipients and controls = NTC cells (ADR-003), clean-control
restriction, D_max 60 um with bins [0, 15], (15, 30], (30, 60] um, HC1 or cluster-robust SEs.
Datasets: Perturb-FISH tumour (all cells) and Perturb-Multi (sections 4, 5, 8, 9, 10).
Simulator gate (run first, results/bench_spatial): E2 with the basis must keep NTC
pseudo-target false positives at p < 0.05 within [0.03, 0.08] in the clonal scenario at every
ring where E2 without the basis or E1 exceed 0.10; otherwise it is not used on real data.
H1' (calibration): NTC pseudo-targets at p < 0.05 in [0.03, 0.08] overall and per ring; at
q < 0.10 at most 0.10.
H2' (autonomous): at least 20% of targets with one autonomous hit at q < 0.10 pooled.
H3' (spillover beyond artifacts): spillover hits at q < 0.10 exceed 3x the NTC-scaled count.
H4' (bleed-through): artifact detector on the E2 autonomous and ring profiles; first-ring R2
minus last-ring R2 is the artifact fraction, reported per dataset whether or not H3' holds.
The real-data runs are launched in parallel with the gate for time reasons; their results are
not read before the gate result is recorded here.

## 2026-10-05 07:43 Audits of Perturb-DBiT (6 samples) and Spatial Perturb-seq (3 chips)

Perturb-DBiT (reports/audit/perturb_dbit): 22,421 pixels, 21.6% with a dominant guide, 526
targets, only 23 NTC pixels over 14 NTC guides, so no NTC calibration is possible and the
pixel-level analysis can only be exploratory. Same-target adjacency z = 25 to 33: pixels of
one tumour clone share the guide, as expected for transduced HT29 or E0771 cells injected and
grown as metastases.

Spatial Perturb-seq (reports/audit/spatial_perturbseq): 229,775 cells, 2.1% with a guide,
17 targets, 112 cells with the single safe-harbour control guide (mSafe), same-target
adjacency z = 28 to 42 (AAV spread from the injection site; see the inspection report). One
control guide with 112 cells across 3 chips allows a weak calibration check only.

Both datasets therefore enter the cross-technology benchmark for autonomous effects and for
the artifact detector, with spillover estimates labelled exploratory where NTC calibration
cannot be run. Clonality analysis for these two did not run (see audit2.log); rerun pending.

## 2026-10-05 07:47 Gate result for E2 with spatial basis (results/a78cc33a76)

30,000 cells on Perturb-FISH geometry, NTC recipients, clean controls, 15% assigned (30% of
them NTC), 2 reps. NTC pseudo-target false positives at p < 0.05, spillover:
- clonal scenario: E2 0.07 / 0.08 / 0.06 and E2+basis 0.07 / 0.07 / 0.06 for rings 0 / 1 / 2;
  E1 identifies nothing in rings 0 and 1 (no clean controls under clonal assignment) and
  0.05 in ring 2.
- spill scenario: ring 0 has only 120 identified tests (one target) and gives 0.15 (E2) and
  0.10 (E2+basis); rings 1 and 2: 0.16 / 0.12 and 0.09 / 0.08.
- no_effect: ring 0 (120 tests) 0.17 / 0.18; ring 1 0.13 / 0.12; ring 2 0.09 / 0.08.
Autonomous: all three estimators at 0.05 to 0.08 with 92 to 95% NTC coverage.
The basis lowers the false discovery proportion where clonal and artifact structure is present
(ring 2 "everything": E2 0.54 vs E2+basis 0.31; autonomous "everything": 0.24 vs 0.11;
clonal autonomous 0.20 vs 0.12) at little cost in power (0.52 vs 0.44, 0.57 vs 0.51).
Gate verdict: passed for rings 1 and 2 under clonal assignment (the criterion's trigger
condition, E1 or plain E2 above 0.10, did not occur there), with the explicit caveat that
ring 0 and the sparse-exposure settings show 0.10 to 0.18 on too few tests (120 to 1,080) to
be conclusive. Rule for the real-data reading: ring-0 results from E2 are reported with this
caveat and are not used for claims; rings 1 and 2 are.

## 2026-10-05 07:47 E2 on real data (pre-registered 07:42)

Perturb-FISH tumour, E2 without spatial basis (results/b78f984897): NTC pseudo-targets at
p < 0.05: 0.162 (autonomous), 0.134 / 0.089 / 0.197 (rings 0 / 1 / 2); at q < 0.10: 0.11,
0.10 / 0.05 / 0.17. Fails H1' everywhere. Same niche confounding as E1 with global strata.

Perturb-FISH tumour, E2 with 40-centre spatial basis (results/246b4deb3e):
- H1': NTC at p < 0.05: 0.058 (autonomous), 0.116 (ring 0), 0.064 (ring 1), 0.070 (ring 2);
  at q < 0.10: 0.017, 0.037, 0.019, 0.016. Passes for autonomous and rings 1 and 2; ring 0
  fails as the gate anticipated (sparse exposure: median 69 NTC recipients per pseudo-target
  and 7.5 per real target).
- H2': 33 of 33 targets have at least one autonomous hit at q < 0.10 (910 hits over 16,500
  tests, 5.5% versus 1.7% for NTC pseudo-targets). Passes. Empirical false discovery
  proportion estimated from the NTC rate: 0.017 / 0.055 = 0.31, so about a third of the
  autonomous hits are expected to be false even after BH, which means the BH q-values are
  not nominal here; the manuscript reports this NTC-estimated FDP next to every hit count.
- H3': spillover hits at q < 0.10: ring 1 837 of 12,000 (7.0%) versus NTC 1.9%, ratio 3.7;
  ring 2 730 of 16,000 (4.6%) versus NTC 1.6%, ratio 2.9; ring 0 363 of 5,000 (7.3%) versus
  NTC 3.7%, ratio 2.0 (not used). Pre-registered threshold 3x: passes at ring 1 (15 to 30 um),
  borderline at ring 2 (30 to 60 um). Estimated FDP among ring-1 spillover hits: 0.27.
- H4' (bleed-through) and the biology follow in the next entry (scripts/analyze_real.py).
Caveat carried from the gate: ring-0 calibration is not established; ring-0 hits are listed
but not claimed.

Perturb-Multi, E2 without basis (results/a2531ecb83): NTC at p < 0.05: 0.216 (autonomous),
0.122 / 0.134 / 0.145; q < 0.10: 0.022, 0.009 / 0.018 / 0.014. Fails H1'. Spillover hits at
q < 0.10 (250 / 199 / 281) do not exceed the NTC-scaled expectation (about 240 / 700 / 500):
no spillover signal. Autonomous: 333 hits in 134 of 201 targets against an NTC rate that
predicts 924 false calls, so the autonomous hit count is below its own null expectation;
the NTC pseudo-target autonomous test (one NTC guide's cells versus the other NTC guides')
is anticonservative here, consistent with guide-specific clonal clusters sitting in their own
niches (same-guide adjacency 260x). The spatial-basis run is in progress.

## 2026-10-05 07:50 H4', biology and prediction on Perturb-FISH (results/246b4deb3e); Perturb-Multi with basis

H4' (bleed-through detector, scripts/analyze_real.py): projection of spillover profiles on the
autonomous profile gives alpha_hat with median -0.04 (ring 0, 13 targets), -0.02 (ring 1, 27),
0.00 (ring 2, 35) and R2 means 0.015 / 0.027 / 0.041. First-ring minus last-ring R2 is -0.04:
no bleed-through signature. The Perturb-FISH spillover profiles are not diluted copies of the
autonomous profiles, so segmentation spill is not what drives the ring-1 excess. Caveat: the
detector only sees proportional transfer; a non-proportional artifact would pass.

Ligand-receptor enrichment: only 8 of 28,000 (target, outcome) pairs are known LR axes
(OmniPath ligrecextra + CellPhoneDB) because the library targets intracellular TLR/NF-kB
signalling and the 500-gene panel is immuno-oncology; 0 of 1,567 hits fall on them.
Uninformative, reported as such.

Niche ranking (rings 1 and 2, score = mean z^2 - 1): MYD88 2.10, MAP3K7 1.91, LBP 1.51,
RELA 1.44, LY96 1.43, TBK1 1.43, JUN 1.43, IKBKB 1.31, TAB2 1.18, TRAM1 1.17; the 95th
percentile of NTC pseudo-target scores is 0.56, so 14 of 33 targets exceed it. The ranking is
dominated by core TLR4/NF-kB adaptors, which is the expected biology of a tumour-cell NF-kB
knockout changing cytokine output to neighbours, and also what a shared niche confound would
produce if the basis under-corrects; the two are separated only by the NTC calibration above.

Held-out-target prediction (kernel ridge, leave-one-target-out): Replogle K562 gwps embedding
covers 25 of 35 targets, GO covers 33. Mean r against observed z-profiles: autonomous 0.09
(null 0.10), spillover 0.21 (null 0.21) for Replogle; 0.07 (0.09) and 0.19 (0.20) for GO.
No embedding predicts either autonomous or spillover profiles beyond a permuted-embedding
null (z between -1.4 and 0.2). Reported as a negative finding (C5): with 33 targets from one
pathway there is nothing to generalise from.

Perturb-Multi, E2 with 40-centre basis (results/6397283d7d): autonomous NTC p < 0.05 fraction
0.072 (passes marginally), spillover 0.135 / 0.128 / 0.129 (fails), q < 0.10: 0.016 and 0.044 /
0.039 / 0.035. Spillover hits (667 / 569 / 605) are below the NTC-scaled expectation (about
1,070 / 1,500 / 1,270); autonomous hits 200 in 103 of 201 targets against an expectation of
670 false calls. No signal beyond the null in Perturb-Multi with this configuration. Forty
centres per 15 mm liver section (about 2 mm spacing) cannot follow lobular zonation at the
25 um scale of hepatocyte neighbourhoods.

Amendment A4 (declared before running): basis-size sensitivity on Perturb-Multi with 150 and
400 centres per section; H1' is re-evaluated per size. Exploratory E2 runs (labelled as such,
no NTC calibration possible) on Perturb-DBiT (unperturbed-pixel controls) and Spatial
Perturb-seq (any_other_guide controls, 112 mSafe cells) for the cross-technology table.

## 2026-10-05 07:56 Amendment A5: E2 ring coefficients from recipients only (bug fix)

The MYD88 case study showed ring-1 effects of -0.04 to -0.12 with q-values below 1e-4 from
45 recipients, which is impossible for a 45-cell mean. Cause: in E2 the ring count m_i^b was
non-zero for cells that carry the target themselves (their clone-mates are their neighbours),
so the "spillover" coefficient was partly estimated from within-clone homogeneity of target
cells, with correspondingly tiny residuals. Fix: m_i^b is set to zero for own-target cells, so
ring coefficients are identified from eligible recipients only (test added). All E2 results
above (results/b78f984897, 246b4deb3e, a2531ecb83, 6397283d7d) are superseded; H1' to H4' are
re-evaluated on the reruns, and the gate benchmark is rerun with the corrected E2. The
superseded directories stay in results/ for the record.

## 2026-10-05 07:56 Clonality in Spatial Perturb-seq and Perturb-DBiT (exploratory)

Spatial Perturb-seq (chip B03018A2 and others; reports/audit/spatial_perturbseq/clonality.md):
all guide-carrying cells are concentrated at the AAV injection site, so pairs of assigned
cells are 5.5 to 6.7 times more frequent than the stratified null for *different* targets
and 8.7 / 7.6 / 6.2 / 5.5 times for the same guide (bins to 15, 30, 50, 75 um). The
guide-specific excess over the shared delivery excess is therefore only 1.3x within 15 um,
decaying to none by 50 um; the misassignment signature is weak (median r 0.17, 1 of 17 above
0.3). Here the dominant structure is delivery, not clonality: spatial strata or a spatial
basis are needed because exposure is confined to one region, but the undetected-sibling
problem is milder than in the tumour screens.

Perturb-DBiT (pixel units; first run used micrometre bins by mistake and is superseded by
reports/audit/clonality5.log): 47.5% of guide-assigned pixels sit in same-guide components,
and unassigned neighbour pixels mirror the perturbed pixels' profile almost perfectly
(median r 0.99 over 20 targets), as expected when a lesion of one clone spans many pixels and
the dominant-guide call fails in some of them. Pixel-level spillover in DBiT is thus mostly
intra-lesion mixing; it is reported as exploratory only.

Correction to the A5 entry: results directories are keyed by the configuration hash only, so
the E2 reruns overwrote results/246b4deb3e, b78f984897, a2531ecb83 and 6397283d7d in place;
the superseded numbers survive only in this notebook (entries of 07:47 and 07:50) and in the
git history of the calibration JSON files. From this commit on, every results directory also
records the git revision that produced it (code_version.txt).

## 2026-10-05 08:03 Corrected E2 (A5): gate fails; amendment A6 (spatially clustered SEs)

Gate rerun with the A5 code (results/a78cc33a76, overwritten in place): spillover NTC false
positives at p < 0.05 are now 0.10 to 0.26 for E2 and E2+basis in rings 1 and 2 under the
clonal and "everything" scenarios (ring 2: 0.24 / 0.22 clonal, 0.26 / 0.23 everything) and
0.11 to 0.16 even in the plain spill scenario; autonomous stays at 0.05 to 0.06; E1 stays at
0.05 to 0.08 where identified. Identifying the ring coefficient from recipients only removed
the spurious precision from clone-mates but exposed the next problem: HC1 standard errors treat
the few recipients around one clone as independent, while their expression is correlated
through the shared niche. The gate criterion (E2+basis within [0.03, 0.08] wherever plain E2
or E1 exceed 0.10 under clonal assignment) is therefore not met, and under the 07:42
pre-registration E2 spillover results are not used for claims.

Corrected Perturb-FISH E2+basis (results/246b4deb3e, A5 code): NTC p < 0.05: autonomous
0.078, rings 0.171 / 0.079 / 0.067; at q < 0.10 spillover hits 2,249 with an NTC-estimated FDP
of 0.53 (pooled over rings; per-ring table printed in the log of this entry). H3' fails.

Amendment A6 (declared before running): standard errors clustered by sample x 100 um spatial
tile (cluster-robust with hundreds of clusters) in E2, point estimates unchanged. Gate rerun
with E1, E2+basis and E2+basis+clustered SEs; the same criterion applies. If A6 passes the
gate, the real-data E2 runs are repeated with it; if not, the manuscript reports E1 with local
tiles as the only calibrated spillover estimator and its null result.

## 2026-10-05 08:06 A6 fails the gate; amendment A7 (permutation inference for E2)

Gate with A6 (results/bench_spatial3.log): spatially clustered SEs raise, not lower, the NTC
false-positive rate of E2+basis in the clonal scenario (ring 2: 0.29 vs 0.22; ring 1: 0.18 vs
0.12); FDP 0.75 to 0.88. Clusters of 100 um around sparse exposed recipients contain too few
cells for the cluster-robust estimator, and the 40-centre basis can absorb part of the exposed
clusters (a basis function centred on a clone fits its recipients), shrinking residuals. Both
analytic SE routes are rejected. E1 (studentized stratified permutation, clean controls) is
calibrated at 0.05 to 0.08 on the same simulated data wherever it is identified.

Amendment A7 (declared before running): E2 gains the same stratified permutation null as E1.
Guide labels are permuted within strata (sample x cell type), exposures recomputed, the
regression refitted per target, and the studentized coefficient t = estimate / HC1 SE is
compared with its permutation distribution; p is the permutation-calibrated z (two-sided
normal on (t_obs - mean t_null) / sd t_null), identification requires at least
max(20, n_perm / 2) valid permutations, and the CI is the HC1 SE rescaled by sd t_null.
Covariates and the spatial basis are kept (they still adjust the point estimate). The gate
criterion is unchanged; if A7 passes, the real-data E2 runs are repeated with it.

## 2026-10-05 08:12 Correction of entry times

The clock times in today's headings from 07:42 onward were first written as estimates and
were wrong (they ran ahead of the clock by up to four hours). They have been replaced by the
times of the git commits that first contained each entry (git log, Europe/Paris). The order
of entries and their content are unchanged.

## 2026-10-05 08:20 A7 gate result (results/b4c9dbe669) and the final configuration

E2 with the stratified permutation null (with or without the spatial basis), 30,000 cells, NTC
recipients, clean controls, 100 permutations, 2 reps:
- autonomous: NTC false positives at p < 0.05 0.04 to 0.06 in every scenario (passes); power
  0.75 to 0.82, FDP 0.09 to 0.20 (E1: 0.12 where identified).
- spillover ring 2: 0.07 (no effect, spill), 0.14 to 0.15 (clonal, everything); E1 0.05 to
  0.08. Ring 1: 0.00 to 0.01 (over-conservative; the studentized null of sparse exposures is
  heavy-tailed relative to the observed statistic), E1 not identified. Power at ring 2 under
  clonal assignment 0.40 to 0.45 against E1's 0.12, at FDP 0.66 to 0.67.
Gate criterion (within [0.03, 0.08] at every ring where plain E2 or E1 exceed 0.10 under
clonal assignment): fails at ring 2 (0.14 to 0.15). Per the 07:42 pre-registration, E2
spillover estimates are not used for claims in any variant (A5, A6, A7). E1 with local tile
strata remains the only calibrated spillover estimator, and on Perturb-FISH it identifies
2,000 to 7,500 tests at 150 to 400 um tiles with zero autonomous hits and 0 to 6 spillover hits.

Final configuration for the manuscript: calibrated inference = E1 with 250 um tiles (primary,
pre-registered tile size) and E1 global as the anticonservative reference; exploratory =
E2 with spatial basis and permutation null, reported with the NTC-scaled expectation and the
implied false discovery proportion next to every count, labelled exploratory. The E2
permutation run on Perturb-FISH follows (configs/e2_perturb_fish_sb40_perm.yaml). Perturb-Multi
E2 permutation is not run: 250 targets on 750,000 cells would take many hours for an
estimator the gate rejected.

## 2026-10-05 08:50 Exploratory E2 with permutation null on Perturb-FISH (results/197e22efb1)

E2, 40-centre basis, NTC recipients, clean controls, 100 stratified permutations:
- NTC pseudo-targets at p < 0.05: 0.053 (autonomous), 0.039 / 0.038 / 0.054 (rings 0 / 1 / 2);
  at q < 0.10: 0.001 and 0.000. On real data the permutation null is calibrated at every
  ring, including the first, which the simulator's clonal scenario did not predict (0.14 at
  ring 2 there); the simulator's clonal structure is harsher than Perturb-FISH's.
- Hits at q < 0.10: 16 autonomous (NTC-scaled expectation 8.2, ratio 1.9) and 2 / 2 / 3
  spillover in rings 0 / 1 / 2 (expectation 0). With a calibrated null, the thousands of E2
  calls of the analytic runs collapse to a handful: the earlier counts were almost entirely
  products of anticonservative standard errors, as the NTC-estimated FDPs had already said.
- Artifact detector on these profiles: first-ring R2 0.030 vs last-ring 0.025 (difference
  0.005), alpha_ring0 0.39 (all 34 targets positive) vs 0.10 in the last ring. A positive,
  first-ring-concentrated projection coefficient with negligible variance explained means a
  small proportional component exists but accounts for about 3% of first-ring spillover
  variance: bleed-through is detectable but minor.
- Ranking: MYD88 remains first (125 cells, 45 NTC ring-1 recipients); scores are now on the
  calibrated scale and are reported with the NTC 95th percentile in the figure.

Final real-data statement: with calibrated inference (E1 tiles or E2 permutation), Perturb-FISH
shows a weak autonomous signal (16 calls against 8 expected) and no spillover signal beyond
the non-targeting expectation at the current sample size; Perturb-Multi shows neither. The
manuscript, figures and tables are regenerated from results/197e22efb1 (Perturb-FISH) and
results/6397283d7d (Perturb-Multi, analytic E2 with basis, exploratory; the 150-centre run is
still in progress and is a sensitivity analysis only).

## 2026-10-05 09:05 Perturb-Multi basis-size sensitivity (amendment A4) and corrected sb40

Both runs use the A5 code (recipients-only ring identification), analytic SEs, NTC recipients,
clean controls. With clean controls only ring 2 (30 to 60 um) reaches positivity in Perturb-Multi.
- 40 centres per section (results/6397283d7d, overwritten by the A5 rerun): 70,015 identified
  tests; NTC p < 0.05: 0.063 (autonomous), 0.096 (ring 2); 223 autonomous and 667 spillover
  calls at q < 0.10 against NTC-scaled expectations that exceed them (FDP 1.0 and 0.62).
- 150 centres (results/fddb2ec60c): NTC p < 0.05: 0.054 (autonomous), 0.032 (ring 2); H1'
  passes for both. 134 autonomous calls against an expectation of about 670 and 154 spillover
  calls against about 180: nothing beyond the null.
The 150-centre run is the Perturb-Multi configuration used for the manuscript's exploratory
numbers (it is the one that passes H1'); the conclusion is identical for both: no autonomous
or spillover signal exceeds the non-targeting expectation in the 209-gene liver panel with
NTC recipients.

## 2026-10-05 08:56 Exploratory E2 on Spatial Perturb-seq (results/df5156edac)

Chip B03018A2 only, 2,000 most detected genes, any_other_guide recipients (one safe-harbour
guide, 112 cells), 40-centre basis, analytic SEs, no clean-control restriction. NTC
pseudo-target rates at p < 0.05: 0.068 (autonomous), 0.148 / 0.054 / 0.059 (rings 0 / 1 / 2).
305 autonomous and 3,014 spillover calls at q < 0.10, both below their NTC-scaled
expectations (FDP 1.0). With 3.4% guide-positive cells, 57% of them at a single guide UMI, and
one control guide, this dataset cannot support a spillover claim; it enters the
cross-technology table as "no signal beyond the null, exploratory".

## 2026-10-05 09:35 Second pass: plan and pre-registration

Framing decided: a methods paper (non-targeting calibration for spatial screens, why
neighbour-effect designs fail, and what sample sizes would be needed), with the biological
null result as the empirical finding. Four work items, in order: (1) robustness and
ablations, (2) power analysis, (3) an undetected-sibling model, (4) release tidying.

Pre-registration, robustness (real data): E1, 250 um tiles, NTC recipients, Perturb-FISH,
with distance bins [0, 10, 20, 40, 60] and [0, 20, 40, 60] in place of [0, 15, 30, 60].
Criterion: the conclusion (NTC p < 0.05 fraction within [0.02, 0.08]; zero autonomous hits;
fewer than 10 spillover hits) is unchanged under both binnings.

Pre-registration, ablation (simulator, Perturb-FISH geometry, 30,000 cells, NTC recipients,
clonal scenario and plain spill scenario, 2 reps): E1 with strata sample x cell type (full),
E1 with strata sample only (no cell type), E1 with 250 um tile strata, E2 with and without
covariates, E0 pseudobulk (no space). Metrics: NTC false positives, power at q < 0.10, FDP.
Expectation written down before running: removing cell type raises NTC false positives where
cell-type composition varies spatially; tiles lower power and keep calibration; removing
covariates changes little in the plain scenario and raises FDP under density effects.

Pre-registration, power analysis (simulator, Perturb-FISH geometry, 30,000 cells, 1 rep per
cell): grid over fraction of cells assigned {0.05, 0.15, 0.30}, fraction of assigned cells
that are NTC {0.1, 0.3, 0.6}, planted spillover LFC sd {0.5, 1.0, 2.0} on 20% of genes for
half of the targets; estimators E1 with 250 um tiles and E2 with spatial basis and permutation
null. Output: power at q < 0.10 in rings 1 and 2 against the median number of NTC recipients
per target, and NTC false positives. Purpose: a design recommendation (how many confirmed
non-target recipients per target a screen needs to detect spillover of a given size), not a
claim about any dataset. No real-data hypothesis is attached.

Undetected-sibling model (exploratory, then tested on the simulator before any real-data use):
the profile of unassigned neighbours of g cells is modelled as a mixture
pi * (autonomous profile of g) + (1 - pi) * (true recipient profile); pi is estimated per
target by projection on the autonomous profile across genes, and the corrected spillover
profile is the residual rescaled by 1 / (1 - pi). On the simulator's clonal and misassignment
scenarios the test is whether the corrected unassigned-recipient estimates have smaller bias
than the uncorrected ones and whether pi tracks the planted misassignment rate.

## 2026-10-05 09:50 Clean-clone check

A fresh `git clone` into a scratch directory followed by `uv sync --extra torch`, `make test`,
`make smoke` and `uv run mypy` succeeded (59 tests passed, 6 skipped for absent raw data,
smoke pipeline produced results/2018dbfee2, mypy clean). The full `make reproduce` was not run
from the clone because it needs about 20 GB of downloads and many hours; its stages are the
same Makefile targets that produced every result in this repository.
