# SpatialSpill pipeline. `make reproduce` regenerates everything from a clean clone.
# Every target is a script in scripts/; outputs land in results/<config_hash>/ or reports/.

UV ?= uv run
PY := $(UV) python
CONFIG ?= configs/default.yaml

.PHONY: help setup lint test smoke data audit reproduce clean-results

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

reproduce: setup data audit ## regenerate every result
	$(PY) scripts/run_pipeline.py --config $(CONFIG)

clean-results:
	rm -rf results/*/
