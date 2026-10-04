# STATUS

Last updated: 2026-10-05 01:30 (paused; background runs continue in tmux)

## Achieved

- Repository, locked environment, CI (lint, mypy, pytest with coverage gate, smoke pipeline
  on a bundled Perturb-DBiT subset), pre-commit hooks, Makefile.
- Dataset access verified and recorded with URLs, sizes and licences (data/README.md).
  Downloaded: Perturb-Multi (14.2 GB), Perturb-DBiT (GSE319277 + Xenium companion files),
  Perturb-map (GSE193460), Perturb-FISH tumour tables (Brain Image Library), Spatial
  Perturb-seq (GSE274447). Five loaders into one schema, each validated.
- Audits for Perturb-Multi, Perturb-FISH, Perturb-map, bundled DBiT (reports/audit/): guide
  call rates, targets, NTC counts, density balance, three neighbour graphs, same-target
  adjacency against a stratified permutation null.
- Finding: same-guide cells are spatially clustered in every in vivo dataset (clonal growth or
  local delivery), and unassigned neighbours of perturbed cells carry the perturbed
  signature. Recipient policy changed accordingly (ADR-003).
- docs/estimands.md (C1): estimands, assumptions, identification, assumption tests, artifact
  signatures; all cited DOIs verified by scripts/check_citations.py.
- Estimators E0 (neighbour t-test, pseudobulk DE baselines), E1 (stratified difference in
  means, studentized permutation null, permutation-calibrated z), E2 (regression adjustment
  with robust or clustered SEs), E3 (augmented IPW with Monte Carlo exposure probabilities),
  E4 (GNN with interventional counterfactuals; ranking only, its SEs are anticonservative).
- Bleed-through artifact detector (projection of first-ring profile on autonomous profile).
- Embeddings (Replogle pseudobulk PCA, GO SVD, ESM2), held-out-target prediction (kernel
  ridge, leave-one-target-out, permuted-embedding null), ligand-receptor enrichment and
  niche-gene ranking modules, all unit-tested; not yet run on real estimates.
- Real-data E1 on Perturb-FISH: three documented amendments. Global strata are
  anticonservative at ring 0 (spatial autocorrelation around clones); 250 um tile strata are
  calibrated but identify only 8,000 of 95,000 tests and find no spillover.
- Simulator benchmark: E1 calibrated on NTC pseudo-targets; E4 best for ranking spillover
  (AUROC about 0.9) but uncalibrated; composition-aware truth table.
- Simulator (C2) with planted autonomous, spillover, bleed-through, density, batch, clonal
  and misassignment effects; benchmark scoring (null FPR, coverage, power, AUROC, FDP).
- Pre-registered first E1 analysis of Perturb-FISH with two documented amendments.

## Running when paused (tmux sessions ss_pow3, ss_e1_tiles, ss_e1_multi)

- results/bench_power2.log: 40,000-cell power benchmark, E0 to E3, 7 scenarios x 2 reps.
- results/e1_fish_tile150.log, e1_fish_tile400.log: Perturb-FISH tile sensitivity.
- results/e1_perturb_multi.log then e1_perturb_multi_tile400.log: pre-registered Perturb-Multi E1.
Run `uv run python scripts/aggregate_results.py` after they finish.

## Not yet achieved

- Per-dataset artifact fractions on real data (detector exists; needs the E1/E2 tables with
  autonomous and first-ring estimates for Perturb-Multi, DBiT, Stereo-seq).
- E2/E3/E4 on real data; cross-dataset benchmark table (C4) with all estimators.
- Held-out gene / held-out technology prediction on real estimates (C5); biology and
  validation (C6); packaging, tutorial, manuscript (C7).
- Loaders for Perturb-DBiT Xenium companion; audits for DBiT full and Stereo-seq.
- Convex-hull or mask-based tissue edge distance; undetected-sibling model; E4 uncertainty
  (retraining ensembles or NTC-based calibration).

## Repository visibility

Created private (ADR-001). To publish: `gh repo edit MMatinGerami/spatialspill --visibility public --accept-visibility-change-consequences`.
