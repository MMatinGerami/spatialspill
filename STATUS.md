# STATUS

Last updated: 2026-10-09 00:40

## Achieved

- Repository, locked environment, CI (lint, mypy, pytest with an 80% coverage gate, about 91% covered, smoke
  pipeline on a bundled Perturb-DBiT subset), pre-commit hooks, Makefile with `reproduce`,
  `real`, `bench`, `analyze` and `paper` targets.
- Five datasets downloaded and loaded into one schema (Perturb-Multi with exact count
  recovery, Perturb-DBiT pixels, Perturb-map Visium, Perturb-FISH tumour, Spatial Perturb-seq
  Stereo-seq), with audits and clonality analyses for all five (reports/audit/).
- docs/estimands.md (C1); simulator and benchmark (C2) with composition-aware truth and NTC
  pseudo-targets as exact nulls; estimators E0 to E4; bleed-through detector (C3);
  embeddings, held-out-target and held-out-technology prediction (C5); niche ranking,
  ligand-receptor enrichment and a case-study script (C6); CITATION.cff, Zenodo metadata, API
  doc, executed tutorial notebook, manuscript draft that compiles with all numbers from
  generated macros (C7).
- Pre-registered real-data analyses with seven declared amendments (A1 to A7), each recorded
  with the failure that motivated it and the superseded numbers.

## Main findings so far

- Same-guide cells are spatially clustered in every in vivo screen (clonal growth or local
  delivery); unassigned neighbours of perturbed cells carry the perturbed signature.
  Recipients must carry a confirmed non-target guide (ADR-003).
- Globally stratified permutation nulls and HC1 regression SEs are anticonservative because
  exposed recipients share a clone's niche. Local tile strata restore E1's calibration but
  remove nearly all power; in Perturb-FISH 0 autonomous and 2 spillover calls survive a local
  null (pre-registered H2 and H3 not supported), and in Perturb-Multi no spillover exceeds the
  non-targeting expectation in any calibrated configuration.
- On the simulator, the published neighbour-t-test design has false discovery proportions of
  0.3 to 0.6; E4 ranks planted spillover best (AUROC about 0.9) but its intervals are unusable.
- Gene embeddings do not predict held-out profiles beyond a permuted null; nothing transfers
  across technologies with the available target sets.

## Second pass (2026-10-05 morning)

- Robustness: tile-size (150, 250, 400 um) and bin sensitivity on Perturb-FISH; ablation on the
  simulator (strata, tiles, covariates, space). Conclusions unchanged; details in NOTEBOOK.
- Power analysis (configs/power_grid.yaml): calibrated detection of large outer-ring spillover
  reaches at most 0.54 power with about 350 confirmed non-target recipients per target on
  30,000 cells; Perturb-FISH has medians of 5.5 (first ring) and 10.5 to 32 (outer rings).
  Design recommendation in the manuscript.
- Undetected-sibling mixture model (src/spatialspill/sibling.py): recovers spillover from
  unassigned neighbours on the simulator; not identified on Perturb-FISH (niche read as
  contamination). Limitation recorded.
- Clean-clone test passed (tests, smoke, mypy). Manuscript extended (ablation, power,
  sibling sections) and recompiles with no undefined references.

## Running or next

- Exploratory Perturb-DBiT E2 finished (results/8d9f6e37b2): not calibrated (NTC q < 0.10
  fraction 0.092, NTC-estimated FDP about 0.4), as expected for pixel units; in the manuscript
  and NOTEBOOK.
- Not done, by design or by data: PyPI release and Zenodo deposit (owner's calls); DepMap and Open Targets checks (no calibrated hit set); sibling model with local
  autonomous profiles; mask-based edge distance; E4 uncertainty.

## Repository visibility

Created private (ADR-001); made public on 2026-10-09 after a pre-publication audit of every
number in the README, STATUS and manuscript against the generated results (corrections in
the NOTEBOOK entry of that date).
