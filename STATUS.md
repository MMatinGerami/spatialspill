# STATUS

Last updated: 2026-10-05

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
- Estimators E1 (stratified difference in means, studentized permutation null), E2
  (regression adjustment with robust or clustered SEs), E3 (augmented IPW with Monte Carlo
  exposure probabilities). E4 (GNN) in progress.
- Simulator (C2) with planted autonomous, spillover, bleed-through, density, batch, clonal
  and misassignment effects; benchmark scoring (null FPR, coverage, power, AUROC, FDP).
- Pre-registered first E1 analysis of Perturb-FISH with two documented amendments.

## Not yet achieved

- Artifact detector (C3) and per-dataset artifact fractions.
- Benchmark table across datasets and estimators (C4); E4; space-blind baselines.
- Prediction from embeddings (C5), biology and validation (C6), packaging and manuscript (C7).
- Convex-hull or mask-based tissue edge distance; undetected-sibling model.

## Repository visibility

Created private (ADR-001). To publish: `gh repo edit MMatinGerami/spatialspill --visibility public --accept-visibility-change-consequences`.
