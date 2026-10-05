# STATUS

Last updated: 2026-10-05 08:36

## Achieved

- Repository, locked environment, CI (lint, mypy, pytest with coverage gate about 90%, smoke
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
  remove nearly all power; in Perturb-FISH no autonomous or spillover effect survives a local
  null, and in Perturb-Multi no spillover exceeds the non-targeting expectation.
- On the simulator, the published neighbour-t-test design has false discovery proportions of
  0.3 to 0.6; E4 ranks planted spillover best (AUROC about 0.9) but its intervals are unusable.
- Gene embeddings do not predict held-out profiles beyond a permuted null; nothing transfers
  across technologies with the available target sets.

## Running or next

- Amendment A7 (permutation inference for E2) implemented and gated: autonomous calibration
  nominal, outer spillover ring 0.14 to 0.15 under clonal assignment (ceiling 0.08): E2 stays
  exploratory. Exploratory E2 permutation run on Perturb-FISH and the 150-centre Perturb-Multi
  run are in tmux (ss_e2perm, ss_e2b); figures, tables and the manuscript are regenerated from
  them with `make paper` when they finish.
- Not done, by design or by data: PyPI release, Zenodo deposit and public visibility (owner's
  calls); DepMap and Open Targets plausibility checks (no calibrated spillover hit set exists to
  check); systematic ablation table; undetected-sibling model; mask-based edge distance;
  E4 uncertainty (ensembles or NTC calibration); Perturb-DBiT Xenium companion loader.

## Repository visibility

Created private (ADR-001). To publish: `gh repo edit MMatinGerami/spatialspill --visibility public --accept-visibility-change-consequences`.
