# SpatialSpill pipeline. `make reproduce` regenerates everything from a clean clone.
# Every target is a script in scripts/; outputs land in results/<config_hash>/ or reports/.

UV ?= uv run
PY := $(UV) python
CONFIG ?= configs/default.yaml

.PHONY: help setup lint test smoke data audit real bench analyze paper reproduce clean-results

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  %-16s %s\n", $$1, $$2}'

setup: ## install the locked environment
	uv sync --extra torch --extra pert

lint: ## ruff + mypy
	$(UV) ruff check . && $(UV) ruff format --check . && $(UV) mypy

test: ## pytest with coverage gate
	$(UV) pytest

smoke: ## pipeline on the bundled subset (what CI runs)
	$(PY) scripts/run_pipeline.py --config configs/smoke.yaml

data: ## download raw data (idempotent, md5-checked)
	$(PY) scripts/download_data.py --all

audit: data ## per-dataset audit reports
	$(PY) scripts/audit_datasets.py --config $(CONFIG)

# Real-data runs reported in the manuscript (final configuration, NOTEBOOK 2026-10-05).
# Pinned so that `make paper` works from a clean clone, where the large estimate files
# are not present; override with FISH_RUN=... after a new `make real`.
FISH_RUN ?= results/197e22efb1
MULTI_RUN ?= results/fddb2ec60c

# Simulator runs shown in the manuscript (NOTEBOOK 2026-10-05). Pinned: picking the newest
# benchmark by date picked up the later ablation and power runs instead.
BENCH_RUN ?= results/a0b119699d
GATE_RUN ?= results/b4c9dbe669

real: ## pre-registered real-data runs (E1 and E2 configs)
	for c in e1_perturb_fish e1_perturb_fish_tile250 e2_perturb_fish_sb0 e2_perturb_fish_sb40 e1_perturb_multi e2_perturb_multi_sb40 e2_perturb_multi_sb150 e2_spatial_perturbseq e2_perturb_dbit; do \
	  $(PY) scripts/run_pipeline.py --config configs/$$c.yaml || exit 1; done

bench: ## simulator benchmarks
	$(PY) scripts/benchmark_simulator.py --config configs/benchmark_small.yaml
	$(PY) scripts/benchmark_simulator.py --config configs/benchmark_power.yaml
	$(PY) scripts/benchmark_simulator.py --config configs/benchmark_spatial.yaml

analyze: ## downstream analysis on the latest E2 runs
	$(PY) scripts/clonality_analysis.py --dataset perturb_fish perturb_map perturb_multi spatial_perturbseq perturb_dbit
	$(PY) scripts/analyze_real.py --results $(FISH_RUN)
	$(PY) scripts/analyze_real.py --results $(MULTI_RUN) --organism mouse
	$(PY) scripts/heldout_technology.py --train $(FISH_RUN) --test $(MULTI_RUN)
	$(PY) scripts/case_study.py --results $(FISH_RUN)
	$(PY) scripts/aggregate_results.py

paper: ## figures, tables and the PDF
	$(PY) scripts/make_figures.py --fish $(FISH_RUN) --multi $(MULTI_RUN) --bench $(BENCH_RUN) --gate $(GATE_RUN)
	cp $(FISH_RUN)/case_top.png paper/figures/fig5_case.png
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex >/dev/null 2>&1 || pdflatex -interaction=nonstopmode main.tex >/dev/null
	$(PY) scripts/check_citations.py paper docs

reproduce: setup data audit real bench analyze paper ## regenerate every result from a clean clone

clean-results:
	rm -rf results/*/
