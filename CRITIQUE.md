# Red-team record

For each phase: the five strongest objections a skeptical reviewer would raise, and for each
whether it was fixed (with evidence) or recorded as a limitation.

Standing checklist (applied every phase): bleed-through, density, batch and tissue edge as
alternative explanations; double-dipping (split samples, pre-registration); multiple testing
across genes x distances x cell types; effect sizes with CIs rather than p-values alone;
robustness to neighbour-graph definition and distance binning; overclaiming causality and
untestable assumptions; undocumented manual steps; a simulator too close to the model's own
assumptions; guide misassignment and low-MOI confounding; comparison fairness to published
tools.

(Entries are added at the end of each phase.)

## Phase 1 (infrastructure and audit), 2026-10-05

Objection 1. "Your unit of analysis is wrong for three of five datasets." Perturb-DBiT is a
microfluidic pixel (20 to 50 um, several cells), Perturb-map is a Visium spot (55 um) with
lesion-level labels, and Stereo-seq cellbin segmentation is imperfect. Resolution: units are
recorded per dataset in `uns["spatialspill"]`, distance bins are set per technology
(scripts/audit_datasets.py SCALES), and no cross-technology table pools estimates across
units; each dataset's estimand is stated with its unit. For Perturb-map the estimand is the
lesion's effect on its microenvironment. Recorded as a limitation: pixel/spot-level
"spillover" cannot separate intra-unit mixing from intercellular effects.

Objection 2. "Same-guide cells cluster (z up to 2,400), so your exposure is not randomised."
Resolution: measured on every dataset (reports/audit/*/clonality.md), interpreted as clonal
growth or local delivery (guide-specific, flat or steeply decaying with distance), and handled
by ADR-003 (recipients must carry a confirmed non-target guide). The simulator has clonal and
misassignment scenarios to show what this does to estimators that ignore it. Still open: an
undetected-sibling model to recover power from unassigned cells.

Objection 3. "Perturb-Multi's deposited matrix is not counts; you may have manufactured
counts." Resolution: the deposited `raw.X` is log1p(normalised) with target_sum 93; the
inverse is exact to 1e-13 for every cell and is unit-tested, with the one documented
ambiguity (a common factor of all counts in a cell is unidentifiable; irrelevant at 209
genes). Counts are used only for size factors and GLMs; nothing depends on the tail of the
distribution beyond what the original study used.

Objection 4. "Row-order linkage in Perturb-FISH is fragile." Resolution: proven, not
assumed: the count table's third column decomposes exactly as a MATLAB linear index of the
coordinates (residual 6e-5 over 187,215 rows), the 500 gene names match 9,432 published rows
one-to-one, and the raw 77-column design agrees with the authors' calls on 98.4% of
single-call cells. The discrepancy (2,547 published cells without a raw guide) is recorded in
reports/audit/perturb_fish_inspection.md; both labels are kept in `obs`.

Objection 5. "Edge distance is a bounding-box proxy and density is a 30 um count; tissue
geometry is richer than that." Resolution: recorded as a limitation. Both covariates enter E2
and E3 only as adjustments; the audit reports density balance between perturbed, NTC and
unassigned cells (balanced in Perturb-Multi: 5.53 / 5.54 / 5.48; Perturb-FISH: 12.98 / 13.02
/ 13.26), which is the quantity that matters for A3. A convex-hull or mask-based edge
distance is listed as future work in STATUS.md.

Objection 6 (added). "Your first permutation p-values could not reject anything, then your
z-calibration produced absurd hits." Resolution: both failures are recorded in NOTEBOOK
(amendments A1, A2) with the mechanism (p floor at 1/201; null sd of a raw difference under
varying exposed-group sizes). The fix is a studentized permutation statistic with a minimum
number of valid permutations, checked on a pure null before use. The discarded hits are not
reported anywhere as findings.

## Phase 2 (estimands and identification), 2026-10-05

Objection 1. "Assumption A4 (no unmeasured spatial confounding) is untestable and, in your
own data, false." Resolution: partly conceded. The Perturb-FISH A2 run showed that globally
stratified permutation nulls are anticonservative because exposed recipients cluster around
one clone and share its niche. Tested implication: NTC pseudo-targets at ring 0 gave 13.4% of
tests at p < 0.05. Mitigation: spatial tile strata (A3) restore calibration (3.6%) at a large
cost in identified tests. The estimand document now states that A4 must be checked with
local nulls, and E2/E3 adjust for density and edge distance but not for unmeasured niche
variables. Limitation recorded.

Objection 2. "Recipients restricted to NTC cells are a tiny, possibly unrepresentative
subset." Conceded as a power limitation (ADR-003). Unassigned cells cannot be recipients
because many are undetected siblings. The `any_other_guide` policy enlarges the pool; its
cost is that recipients carry their own autonomous effects, which strata and covariates do not
remove and which the simulator's any_other_guide runs quantify.

Objection 3. "Your exposure mapping (ring counts of g cells) ignores how many cells of *other*
targets surround the recipient." Resolution: E1 and E3 use the clean-control restriction (no
other perturbed neighbour) where positivity allows; E2 can include other-target counts as
covariates; E4 models all guides jointly. Where the restriction is impossible (Perturb-DBiT,
23% of pixels perturbed) it is dropped and said so.

Objection 4. "The non-targeting calibration is only as good as the NTC guides: control
guides can have effects (cutting, immune response to Cas9 expression)." Resolution: NTC
pseudo-targets are tested against the other NTC guides, so a shared Cas9 or delivery effect
cancels; a guide-specific control effect would appear as a calibration failure for that
guide and is reported per guide (the pipeline table has one row per NTC guide).

Objection 5. "A positivity failure is silently set to zero somewhere." Resolution: every
(target, ring, cell type) with fewer than min_cells recipients or controls, or fewer than
min_valid_perm valid permutations, is flagged `identified = False`, excluded from BH, and
counted in the calibration JSON (`n_tests`). The number of identified tests per run is in
every results directory and in the notebook (95,000 global vs 8,000 at 250 um tiles).

## Phase 3 (estimators and simulator), 2026-10-05

Objection 1. "The simulator shares the estimators' assumptions, so good benchmark numbers are
circular." Resolution: the simulator generates negative-binomial counts with multiplicative
effects, exponential-kernel dose-dependent spillover, transcript relocation for bleed-through,
density and batch log-offsets, clonal assignment and barcode misassignment; the estimators
work on log-normalised means with ring indicators or counts and know nothing about kernels.
The benchmark exposed estimator failures (anticonservative ring 0, compositional bias of
normalisation, E3 anticonservative with small exposed groups, E4 uncalibrated), which is what
a non-circular simulator should do. Remaining gap: no simulated segmentation errors beyond
proportional transfer, no cell-type-specific responders in the default scenarios.

Objection 2. "You changed the test statistic twice after seeing real data." Conceded and
recorded (amendments A1, A2, A3; A4 for E2). Each change was declared before the rerun, the
discarded results are kept in `results/` and described in NOTEBOOK, and each statistic's
calibration was checked on a pure null before use. The honest summary is that the first
pre-registered design (E1, exact permutation p, global strata) was both underpowered and
confounded, and the project documents the path to a calibrated estimator rather than hiding
it.

Objection 3. "Your q-values are not nominal: the NTC pseudo-target hit rate at q < 0.10 is
1.6 to 4.4%, and you estimate FDPs of 0.27 to 0.31 among hits." Conceded. BH assumes p-values
that are uniform under the null; the residual niche confounding leaves them mildly
anticonservative. The manuscript reports, next to every hit count, the NTC-scaled expectation
and the implied FDP, which is the empirical null the data themselves provide. Hits are never
listed without that column.

Objection 4. "The GNN is a black box with made-up error bars." Conceded; E4's intervals are
not reported as inference (NOTEBOOK 2026-10-05, benchmark: null FPR above 0.9). E4 is used for
ranking only, where it outperforms E1 and E2 on the simulator (AUROC about 0.9 versus 0.6 at
10,000 cells). Fixing its uncertainty (retraining ensembles or NTC calibration) is future work.

Objection 5. "Comparison with published tools is unfair: you implemented the neighbour t-test
yourself." Partly conceded. E0 reproduces the design of the published analyses (control cells
with versus without a perturbed neighbour, per-gene test, no strata, no permutation) rather
than their exact code (Perturb-FISH used FR-Perturb on neighbour design matrices). On the
simulator this design has the highest apparent power and FDP 0.27 to 0.59; its NTC false
positive rate reaches 10% under clonal assignment. The manuscript states that E0 is a
re-implementation of the design, not the authors' software.

## Phase 4 (prediction and generalisation), 2026-10-05

Objection 1. "Thirty-three targets from one pathway cannot support a held-out-gene claim
either way." Conceded; this is why the result is reported as a failure to beat a
permuted-embedding null, not as evidence that embeddings are uninformative. The code and the
protocol (leave-one-target-out kernel ridge, permuted-embedding null) are in place for larger
libraries; Perturb-Multi's 200 targets give a second, mouse-liver test with the same outcome.

Objection 2. "Human and mouse symbols are aligned by upper-casing, which is not orthology."
Conceded as a simplification (docs/api.md, embeddings module docstring). For the datasets at
hand it affects the GO embedding only; the Replogle embedding is human-only and covers 25 of
35 Perturb-FISH targets.

Objection 3. "Held-out technology with zero overlapping targets is not a transfer test." It
tests whether a niche-score predictor learned on one screen's targets ranks another screen's
targets through a shared embedding. With no overlap and no signal in the test set it cannot
succeed, and the empirical p of 0.11 against the permuted null says exactly that. The honest
reading is "untestable with public data", recorded as such.

Objection 4. "You used E2 z-profiles as the prediction target while E2 was not calibrated."
Partly conceded. The response matrix uses z-statistics as a continuous summary, so
miscalibration scales but does not reorder them; the null-vs-observed comparison is
unaffected. The runs are repeated with the final estimator configuration.

Objection 5. "No ablations." Partly conceded. Removing space (E0 pseudobulk), cell type
(strata = sample only) and embeddings (permuted null) are each present somewhere in the
pipeline, but not as one systematic ablation table; listed in STATUS as not achieved.

## Phase 5 (biology and validation), 2026-10-05

Objection 1. "The niche ranking recapitulates the library (TLR/NF-kB adaptors), so it is
either trivially right or trivially confounded." Conceded. A tumour-cell NF-kB knockout
changing cytokine output to neighbours is the expected biology, and a clone-shared niche would
produce the same ranking. Only the non-targeting calibration separates them, and that
calibration is not established at the first ring and not met by E2 in the simulator under
clonal assignment (gate failures recorded). The ranking is therefore presented as a
hypothesis list with its NTC-estimated false discovery proportion, not as validated biology.

Objection 2. "Ligand-receptor enrichment with 8 testable pairs is meaningless." Conceded;
reported as uninformative. The library and panel do not contain the ligand-receptor space.

Objection 3. "No independent validation against Perturb-map immune phenotypes." Conceded in
part: Perturb-map's named knockouts (Tgfbr2, Ifngr2, Jak2) do not overlap the Perturb-FISH
library (TLR/NF-kB) or the Perturb-Multi liver library, and the lesion-level unit does not
match the cell-level estimand. DepMap and Open Targets plausibility checks were not run
because no spillover gene list reached a calibrated hit set; listed in STATUS as not
achieved.

Objection 4. "The proposed wet-lab experiment cannot distinguish spillover from niche
either." It can: a co-culture with randomised mixing breaks the clonal spatial confound that
the in vivo data cannot, and the conditioned-medium arm separates secreted from
contact-dependent effects. What it cannot do is reproduce the in vivo niche; that is the
trade-off stated in the manuscript.

Objection 5. "One case study is anecdote." Agreed; it illustrates the output format and the
pre-specified readout for validation. Its target changed between the A2 and A5 runs (MYD88 to
MAP3K7), which is itself evidence of how fragile the ranking is at this sample size.
