# Lab notebook

Dated entries, newest at the bottom. Attempts, failures, pivots and open questions are
recorded here as they happen. Pre-registrations appear under a dated "Pre-registration"
heading before the corresponding real-data analysis is run; anything not pre-registered is
labelled exploratory.

## 2026-10-04 Project start

- Environment: Python 3.11 via uv; squidpy 1.8.2, scanpy 1.11.5, anndata 0.12.19,
  mudata 0.3.10, pertpy 1.0.3, torch 2.14.1 (MPS available), torch_geometric 2.8.0.
- Failure: pertpy 1.0.3 fails to import with statsmodels 0.15.0 (`multipletests` moved).
  Fix: pin statsmodels <0.15 (ADR-002).
- Dataset access verification launched in parallel for Perturb-FISH, Perturb-DBiT
  (GSE319277), Perturb-map, Perturb-Multi, other spatial screens, and supporting resources
  (Replogle via scPerturb, atlases, ligand-receptor databases, embeddings). Results recorded
  in data/README.md when they arrive.
- Compute budget: one Apple Silicon machine, 24 GB unified memory, 15 cores, no CUDA. Heavy
  steps will be subsampled and the downscaling documented.
