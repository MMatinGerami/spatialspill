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
