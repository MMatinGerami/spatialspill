# spatialspill API

Install: `pip install spatialspill` (after release) or `uv sync` from a clone.

## Data model

`spatialspill.schema`: one AnnData per dataset with raw counts in `X`, coordinates in
`obsm["spatial"]` (micrometres unless `uns["spatialspill"]["units"]` says otherwise) and the
obs columns `guide`, `guide_confidence`, `target`, `is_ntc`, `is_perturbed`, `area`,
`cell_type`, `sample`, `batch`, `edge_distance`. `validate(adata)` checks it.

Loaders (`spatialspill.loaders`): `load_perturb_multi`, `load_perturb_fish_tumor`,
`load_perturb_dbit`, `load_perturb_map`, `load_spatial_perturbseq`, dispatch through
`loaders.registry.load_dataset(name, **kwargs)`. Each needs the raw files listed in
`data/README.md` (download with `scripts/download_data.py`, `scripts/download_perturb_fish.sh`).

## Exposure and groups

`exposure.compute_exposure(adata, bins_um)` returns ring counts of every target (and NTC
pseudo-target) around every cell within the last bin edge. `estimators.groups.GroupConfig`
sets the recipient policy (`ntc`, `any_other_guide`, `unperturbed`) and the clean-control
restriction; `group_masks` yields the treated, ring-exposed and control indicators.

## Estimators

All share `fit(adata, exposure, outcomes, outcome_names) -> EstimateTable`, whose `.df` has one
row per (target, kind, ring, cell type, outcome) with `estimate`, `se`, `ci_low`, `ci_high`,
`pvalue`, `n_treated`, `n_control`, `identified`; `.with_fdr()` adds BH `qvalue` over
identified rows.

- `NeighbourTTest`, `PseudobulkDE` (E0 baselines)
- `E1Stratified(n_perm, strata_keys, min_cells, groups, pvalue="z"|"perm", ci="perm"|"analytic")`
- `E2GLM(strata_keys, covariates, groups, spatial_basis, cluster_by_sample)`
- `E3DoublyRobust(n_draws, groups, pi_floor, n_boot)`
- `E4GNN(epochs, groups, device)` (ranking only; see its docstring on uncertainty)

## Simulation and scoring

`simulator.synthetic_geometry`, `simulator.geometry_from(adata)`, `simulator.simulate(geom,
SimConfig)`, `simulator.default_scenarios(seed, base)`; `benchmark.score(table, sim_adata)`.

## Artifacts, embeddings, prediction, biology

`artifacts.bleedthrough_table`, `artifacts.artifact_fraction`, `artifacts.covariate_shift`;
`embeddings.perturbseq_embedding`, `go_embedding`, `esm2_embedding`;
`prediction.response_matrix`, `prediction.loto_predict`; `biology.load_lr_pairs`,
`lr_enrichment`, `rank_niche_genes`.

## Scripts

`scripts/run_pipeline.py --config configs/<name>.yaml` (real data, E1 or E2),
`scripts/benchmark_simulator.py --config configs/benchmark_*.yaml`,
`scripts/audit_datasets.py`, `scripts/clonality_analysis.py`, `scripts/aggregate_results.py`,
`scripts/check_citations.py`. Results are keyed by config hash under `results/<hash>/`.
