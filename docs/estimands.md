# Estimands and identification for spillover in spatial CRISPR screens

Written before the estimator code (Phase 2 requirement). Notation is fixed here and reused in
`src/spatialspill/estimands.py`, the manuscript and the tests.

## 1. Setting

A spatial screen yields cells i = 1, ..., n. Cell i has a sample s(i) (slide, field of view
or animal), a cell type c(i), a position x_i in the plane, a segmentation area a_i, a distance
to the tissue edge b_i, and an outcome vector Y_i (expression of G genes, or a derived
phenotype such as a programme score). The treatment is the guide call

    Z_i in {0} U {1, ..., K},

where Z_i = 0 means "no targeting guide" and k indexes the target genes. Cells with a
non-targeting control (NTC) guide are the reference group: they went through delivery and
selection but carry no intended edit. Cells with no guide call are kept as a secondary
reference and never pooled with NTC without saying so.

Let z = (z_1, ..., z_n) be the full assignment vector and Y_i(z) the potential outcome of cell
i under z. Without restrictions there are (K+1)^n potential outcomes per cell and nothing is
identifiable. The restrictions below are the content of this document.

## 2. Interference structure

Assumption A1 (partial interference, Hudgens and Halloran 2008). Cells in different samples do
not interfere: Y_i(z) = Y_i(z_{S(i)}), where z_{S(i)} is the assignment within i's sample.

Assumption A2 (local exposure mapping, Aronow and Samii 2017). There is a maximal radius
D_max and an exposure mapping

    e_i(z) = ( z_i , m_i^1(z), ..., m_i^B(z) ),    m_i^b(z) = sum_{j != i} 1{ d_ij in bin b } 1{ z_j = g }

for a target gene g under study and distance bins b = 1..B partitioning (0, D_max], such that
Y_i(z) = Y_i(e_i(z)). In words: a cell's outcome depends on its own guide and on how many
cells carrying guide g sit in each distance ring around it, and on nothing else about z.

The same mapping with 1{z_j = g} replaced by an indicator 1{m_i^b >= 1} gives the binary
exposure used by the simplest estimator (E1). Cells whose neighbours carry *other* targeting
guides are a nuisance exposure; at the guide multiplicities in real screens (one perturbed
neighbour is already rare for most genes), E1 to E3 restrict recipients to cells with no
perturbed neighbour of any other gene within D_max ("clean controls") and the manuscript
reports how much data that removes. E4 models all guides jointly.

## 3. Estimands

All expectations are over the (hypothetical) re-randomisation of guides within strata and the
population of cells in the stratum.

Cell-autonomous effect of gene g in recipient type c:

    tau_g(c) = E[ Y_i(z_i = g, m_i = 0) - Y_i(z_i = 0, m_i = 0) | c(i) = c ]

Both potential outcomes have no perturbed neighbour within D_max, so tau_g(c) is the effect of
being perturbed in an otherwise unperturbed neighbourhood.

Spillover (non-cell-autonomous) effect of gene g at distance bin b on recipient type c:

    tau_g(b, c) = E[ Y_i(z_i = 0, m_i^b >= 1, m_i^{-b} = 0) - Y_i(z_i = 0, m_i = 0) | c(i) = c ]

The recipient is unperturbed, has at least one g-perturbed cell in ring b and none in other
rings, versus no g-perturbed cell within D_max. The per-neighbour (dose) version replaces the
indicator by a linear term in m_i^b and is what E2 and E3 report by default; E1 reports the
indicator version. Both are in the output tables with their definition.

Aggregate spillover of gene g: tau_g^spill = sum over (b, c) weighted by the stratum sizes of
|tau_g(b, c)| / se, used only for ranking (C6), never for inference.

Outcome scale. Y is log1p of size-factor-normalised counts for per-gene estimands and the
raw count scale (negative binomial GLM) for E2. The scale is in every table.

## 4. Identification

Assumption A3 (ignorable assignment within strata). Within a stratum (sample x cell type) the
guide vector is exchangeable with respect to the potential outcomes:

    Z_{S} independent of { Y_i(.) : i in S } | sample, cell type.

This is what lentiviral delivery at low MOI is designed to deliver. Known violations:
cell-type-dependent infectivity (handled by stratifying on cell type), density-dependent
infectivity (checked: covariate balance of local density across guides), and clonal
expansion after delivery (perturbed cells cluster because they are siblings, not because they
were placed; see Section 6).

Assumption A4 (no unmeasured spatial confounding beyond covariates). Conditional on the
strata and the covariates W_i = (local density, area, edge distance), the exposure
e_i is independent of the potential outcomes. This is the key untestable assumption; what is
testable is its observable implication under the design: the distribution of W across
exposure levels inside a stratum must match what re-randomising guides within the stratum
produces.

Assumption A5 (positivity). Every (stratum, exposure level) that enters an estimand has a
positive probability under the design; empirically, at least n_min cells.

Proposition (identification). Under A1, A2, A3, A4 and A5, for any exposure level e,

    E[ Y_i | e_i = e, stratum, W_i ] = E[ Y_i(e) | stratum, W_i ],

so tau_g(c) and tau_g(b, c) are contrasts of conditional means, estimable by (i) stratified
differences of means when W is balanced by design, (ii) regression adjustment, (iii) inverse
probability weighting with exposure probabilities pi_i(e) = P(e_i(Z) = e) computed exactly by
re-randomising Z within strata on the fixed geometry (Aronow and Samii 2017), or (iv) the
augmented (doubly robust) combination of (ii) and (iii). E1 is (i), E2 is (ii), E3 is (iv).

Positivity failure (A5) is reported per dataset as the fraction of (gene, bin, cell type)
cells that fall below n_min; those estimands are marked "not identified", not set to zero.

## 5. Tests of assumptions

| Assumption | Observable implication | Test in this project |
|---|---|---|
| A1 | none within a dataset | respected by building graphs within samples |
| A2 (D_max) | estimated effects in rings beyond D_max are null | far-field calibration: bins beyond D_max must give the same false-positive rate as NTC |
| A3, A4 | NTC guides have zero autonomous and zero spillover effect | mandatory NTC calibration: each NTC guide is treated as a "gene" and run through every estimator; the empirical false-positive rate at FDR q must be at most q; coverage of CIs for the zero effect must be at least nominal |
| A3, A4 | covariate balance | within strata, standardised mean differences of W between exposure levels compared with their permutation distribution |
| A3 (clonality) | guide labels are spatially autocorrelated within strata | join-count statistic for same-guide adjacency against the stratified permutation null; datasets with clonal structure are flagged and analysed with lesion-level units |
| A5 | counts per exposure cell | reported table of identified vs not identified estimands |

## 6. Artifacts that mimic spillover

Bleed-through (segmentation spillover of transcripts). If segmentation assigns part of
cell i's transcripts to its neighbour j, then the measured outcome is
Y_j^obs = Y_j + alpha_ij Y_i with alpha_ij decaying within about one cell radius. The
apparent spillover of gene g onto adjacent cells is then alpha times the *autonomous* effect
vector of g across all genes. Signature: (a) confined to the first ring, (b) the neighbour
effect vector across genes is proportional to the autonomous effect vector, (c) strongest for
gene g's own transcript when the perturbation changes it (knockdown or NMD). The detector
projects the first-ring spillover profile onto the autonomous profile; the projected part is
the bleed-through estimate, the residual is candidate biology, and the per-dataset artifact
fraction is the share of first-ring spillover variance explained by the projection (C3).

Density. Perturbations that change proliferation or survival change local density, and density
changes neighbours' expression. This is a real indirect effect but not the signalling
interpretation most readers attach to "spillover". All estimators are run with and without
local density in W and the difference is reported.

Tissue edge. Cells near the edge have fewer neighbours, different exposure probabilities and
distinct expression. Edge distance is in W and the audit reports effects as a function of b_i.

Batch. Samples are strata; nothing is estimated across samples; sample-level random effects
are in E2.

Clonal structure. Where perturbed cells form clones (Perturb-map lesions), a cell's neighbours
are mostly its siblings and the design is not a cell-level randomisation. The unit becomes the
lesion and the estimand becomes the effect on the lesion's non-tumour microenvironment, with
lesions as the randomised units.

## 7. Multiple testing and reporting

Tests are indexed by (gene, bin, cell type) within an estimator and dataset; BH-FDR is applied
over that whole index set per dataset, and the number of tests is in each table. Every
estimate is reported with a 95% interval (permutation or sample-level bootstrap as stated
per estimator), never a p-value alone.

## 8. References (DOIs checked by scripts/check_citations.py)

- Hudgens MG, Halloran ME. Toward causal inference with interference. J Am Stat Assoc 2008.
  doi:10.1198/016214508000000292
- Aronow PM, Samii C. Estimating average causal effects under general interference, with
  application to a social network experiment. Ann Appl Stat 2017. doi:10.1214/16-AOAS1005
- Phipson B, Smyth GK. Permutation P-values should never be zero. Stat Appl Genet Mol Biol
  2010. doi:10.2202/1544-6115.1585
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. J R Stat Soc B 1995.
  doi:10.1111/j.2517-6161.1995.tb02031.x
