# Decisions (ADR style)

Each entry: context, decision, alternatives considered, rationale, consequences.

## ADR-001 Repository visibility at creation (2026-10-04)

Context: the project brief asks for `gh repo create spatialspill --public`. Matin's standing
workstation policy says repositories are created private and made public only on his explicit
go in the current conversation, and his other portfolio repositories were kept private until
he reviewed them.

Decision: create the repository private, with the brief's README attribution section included
as written, and record the one-line command to flip visibility in STATUS.md. Matin can make it
public after review.

Alternatives: create public as the brief says.

Rationale: public-to-private is not clean (indexed content persists); private-to-public is one
command. The brief's intent (a public, honest record) is preserved and only delayed to review.

Consequences: CI runs on a private repo (free minutes apply). Nothing else changes.

## ADR-002 Toolchain (2026-10-04)

Decision: Python 3.11 via uv with a lockfile; scverse stack (anndata 0.12, mudata 0.3,
scanpy 1.11, squidpy 1.8, pertpy 1.0); statsmodels pinned below 0.15 because pertpy 1.0.3
imports `multipletests` from a statsmodels sandbox path removed in 0.15; torch 2.14 with MPS
and torch_geometric 2.8 for E4; Makefile pipeline (snakemake not installed; Makefile keeps the
dependency graph readable and `make reproduce` is required by the brief); Hydra/OmegaConf
configs hashed to key results directories.

Alternatives: snakemake (adds a dependency, same outcome); PyMC instead of numpyro (numpyro
is lighter on Apple Silicon).

## ADR-003 Recipients must carry a confirmed non-target guide (2026-10-05)

Context: in every in vivo dataset, same-guide cells are spatially clustered (clonal growth or
local delivery), and in Perturb-FISH unassigned cells adjacent to g cells carry g's autonomous
expression signature (median r = 0.72). With guide-call rates of 4 to 5%, "no guide call" is
uninformative about perturbation status.

Decision: spillover recipients and controls are cells with a confirmed guide different from
the target (default: NTC cells; option: any other guide). Cells without a call are excluded
from recipient and control roles. The docs/estimands.md Section 3 definition is updated.

Alternatives: use all unassigned cells (more power, confounded); model the probability of
being an undetected sibling (future work, needs a detection-efficiency model).

Consequences: power is limited by the number of NTC cells near perturbed cells; positivity
failures are reported explicitly.
