# SpatialSpill

Calibrated estimation of cell-autonomous and non-cell-autonomous (spillover) effects in
spatial CRISPR screens, with artifact detection and cross-technology benchmarking.

Status: in progress. See [STATUS.md](STATUS.md) for what is done and what is not,
[NOTEBOOK.md](NOTEBOOK.md) for the dated lab notebook, [DECISIONS.md](DECISIONS.md) for
design decisions and [CRITIQUE.md](CRITIQUE.md) for the red-team record.

## Question

In a spatial CRISPR screen, knocking out gene g in cell i changes cell i itself
(cell-autonomous effect, tau_g) and may change nearby unperturbed cells (spillover,
tau_g(d, c), as a function of distance d and recipient cell type c). Published screens report
intercellular effects with ad hoc statistics. This project asks: which of those effects
survive a formal treatment of interference, calibration against non-targeting guides, and
explicit modelling of segmentation bleed-through, density, batch and tissue-edge artifacts?

## Contributions (planned)

- C1 Estimands and identification for spillover under network interference
  ([docs/estimands.md](docs/estimands.md)).
- C2 A ground-truth simulator with planted autonomous, spillover, bleed-through, density and
  batch effects; estimators benchmarked on bias, interval coverage, power and AUROC.
- C3 An artifact detector separating true spillover from segmentation bleed-through, density,
  batch and tissue-edge effects, with a per-dataset artifact fraction.
- C4 A benchmark of estimators E1 to E4 plus space-blind baselines across at least three
  spatial-screen technologies.
- C5 Held-out-gene and held-out-technology prediction of spillover from gene embeddings.
- C6 A ranked, independently checked list of niche-remodelling genes with mechanistic
  hypotheses and proposed wet-lab validation.
- C7 An installable `spatialspill` package, tutorial, manuscript draft and Zenodo metadata.

## Layout

```
src/spatialspill/   package (loaders, graphs, estimators, simulator, artifact detector)
scripts/            every figure, table and number is produced here
configs/            Hydra/YAML configs; results are keyed by config hash
results/<hash>/     outputs of scripts (tables, figures, metrics)
reports/audit/      per-dataset audit reports
docs/               estimands, methods notes, tutorial
data/raw            immutable downloads (not tracked; see data/README.md)
data/bundled        tiny subset tracked for CI and the smoke pipeline
paper/              LaTeX manuscript draft
tests/              pytest suite
```

## Reproduce

```
uv sync --extra torch
make reproduce      # downloads data, runs every phase, regenerates results/ and figures
make smoke          # runs the pipeline on the bundled subset (what CI runs)
```

## How this project was built

This project was built end-to-end by Claude Code under the direction of Matin Gerami, who
defined the research question and reviewed the outputs. The git history is a true record of
how the work happened; commit dates and author information were not altered.
