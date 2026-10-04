"""Ground-truth simulator for spatial perturbation screens (C2).

Geometry (coordinates, sample, cell type, area, edge distance) is taken from a real dataset
or generated; guides are assigned by one of several designs; outcomes are negative-binomial
counts with planted effects:

- autonomous: perturbed cell i with target g gets log-fold change ``tau_auto[g, :]`` on its
  own expression;
- spillover: unperturbed cell j gets ``tau_spill[g, :] * sum_i w(d_ij)`` over g-perturbed
  neighbours i, with kernel ``w(d) = exp(-d / spill_scale)`` truncated at ``spill_max``;
- bleed-through: a fraction ``alpha`` of each cell's *expected* transcripts is relocated to
  immediately adjacent cells (shared Delaunay edge within ``bleed_max``), so a neighbour of a
  perturbed cell inherits a diluted copy of the autonomous effect;
- density: a log-linear effect of local density on all genes (``beta_density``);
- batch: per-sample log offsets per gene;
- clonality: guide assignment can be clustered (seed cells plus neighbours within a radius
  take the same guide), reproducing the same-guide adjacency excess observed in real data;
- misassignment: with probability ``p_misassign`` an unassigned neighbour of a perturbed cell is
  *labelled* with that guide without receiving its effect (barcode spill).

Everything is stored in ``adata.uns["truth"]`` so benchmarks can score bias, coverage, power
and AUROC. The simulator deliberately does not share code with the estimators: its outcome
model is a count model with multiplicative effects, while E1 works on log-normalised means and
E2 fits a GLM; agreement is therefore a test, not a tautology.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.spatial import cKDTree

from spatialspill.graphs import GraphConfig, build_graph
from spatialspill.schema import SpatialScreenSchema, edge_distance_from_coords


@dataclass
class SimConfig:
    n_genes: int = 100
    n_targets: int = 20
    n_ntc_guides: int = 10
    frac_assigned: float = 0.08  # fraction of cells with a guide call
    frac_ntc_of_assigned: float = 0.2
    # effects (log scale); a fraction of (target, gene) pairs is non-zero
    frac_auto_nonzero: float = 0.1
    auto_lfc_sd: float = 1.0
    frac_spill_targets: float = 0.3  # fraction of targets with any spillover
    frac_spill_genes: float = 0.05  # genes affected per spilling target
    spill_lfc_sd: float = 0.5
    spill_scale_um: float = 20.0
    spill_max_um: float = 60.0
    spill_cell_types: list[str] | None = None  # recipient types that respond (None = all)
    # artifacts
    alpha_bleed: float = 0.0  # fraction of transcripts relocated to touching neighbours
    bleed_max_um: float = 15.0
    beta_density: float = 0.0  # log-linear density effect per neighbour (within 30 um)
    batch_sd: float = 0.0
    clonal_radius_um: float = 0.0  # 0 = independent assignment
    clonal_size: int = 1
    p_misassign: float = 0.0
    # counts
    mean_counts_per_cell: float = 300.0
    nb_dispersion: float = 0.3  # 1/size; variance = mu + phi mu^2
    seed: int = 0


def synthetic_geometry(
    n_cells: int = 20000,
    n_samples: int = 2,
    field_um: float = 1500.0,
    cell_types: tuple[str, ...] = ("A", "B", "C"),
    seed: int = 0,
) -> AnnData:
    """Poisson-ish geometry with a mild cell-type spatial gradient and a hole (tissue edge)."""
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, field_um, size=(n_cells, 2))
    sample = np.array([f"s{i % n_samples}" for i in range(n_cells)])
    # cell type depends on x to create spatial structure
    p = np.stack(
        [
            np.exp(-((xy[:, 0] / field_um - c) ** 2) / 0.1)
            for c in np.linspace(0.2, 0.8, len(cell_types))
        ],
        1,
    )
    p /= p.sum(1, keepdims=True)
    ct = np.array([rng.choice(cell_types, p=pi) for pi in p])
    obs = pd.DataFrame(
        {
            "cell_type": ct,
            "sample": sample,
            "batch": sample,
            "area": rng.lognormal(np.log(150), 0.3, n_cells),
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    obs["edge_distance"] = edge_distance_from_coords(xy, sample)
    a = AnnData(obs=obs)
    a.obsm["spatial"] = xy
    return a


def geometry_from(adata: AnnData, n_cells: int | None = None, seed: int = 0) -> AnnData:
    """Copy coordinates, sample, cell type, area and edge distance from a real dataset."""
    rng = np.random.default_rng(seed)
    idx = np.arange(adata.n_obs)
    if n_cells is not None and n_cells < adata.n_obs:
        idx = np.sort(rng.choice(idx, n_cells, replace=False))
    sub = adata[idx]
    obs = sub.obs[["cell_type", "sample", "batch", "area", "edge_distance"]].copy()
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["sample"] = obs["sample"].astype(str)
    obs["batch"] = obs["batch"].astype(str)
    a = AnnData(obs=obs)
    a.obsm["spatial"] = np.asarray(sub.obsm["spatial"], dtype=float)
    return a


def _assign_guides(
    geom: AnnData, cfg: SimConfig, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = geom.n_obs
    targets = [f"T{k:03d}" for k in range(cfg.n_targets)]
    ntcs = [f"NTC_{k:02d}" for k in range(cfg.n_ntc_guides)]
    guide_pool = [f"{t}_g{j}" for t in targets for j in (1, 2)] + ntcs
    p_ntc = cfg.frac_ntc_of_assigned
    w = np.array(
        [(1 - p_ntc) / (2 * cfg.n_targets)] * (2 * cfg.n_targets)
        + [p_ntc / cfg.n_ntc_guides] * cfg.n_ntc_guides
    )
    label = np.full(n, "none", dtype=object)
    n_assigned = round(cfg.frac_assigned * n)
    xy = geom.obsm["spatial"]
    sample = geom.obs["sample"].to_numpy()
    if cfg.clonal_radius_um <= 0 or cfg.clonal_size <= 1:
        cells = rng.choice(n, n_assigned, replace=False)
        label[cells] = rng.choice(guide_pool, n_assigned, p=w)
    else:
        n_seeds = max(1, n_assigned // cfg.clonal_size)
        seeds = rng.choice(n, n_seeds, replace=False)
        assigned = 0
        for s_idx in seeds:
            if assigned >= n_assigned:
                break
            g = rng.choice(guide_pool, p=w)
            tree_mask = sample == sample[s_idx]
            cand = np.flatnonzero(tree_mask)
            tree = cKDTree(xy[cand])
            near = cand[tree.query_ball_point(xy[s_idx], cfg.clonal_radius_um)]
            near = near[label[near] == "none"]
            take = near[: cfg.clonal_size] if len(near) >= cfg.clonal_size else near
            label[take] = g
            assigned += len(take)
    is_ntc = np.array([str(g).startswith("NTC_") for g in label])
    target = np.array(
        [
            "NTC" if nt else (str(g).rsplit("_g", 1)[0] if g != "none" else "none")
            for g, nt in zip(label, is_ntc)
        ],
        dtype=object,
    )
    return label, target, is_ntc


def simulate(geom: AnnData, cfg: SimConfig | None = None) -> AnnData:
    cfg = cfg or SimConfig()
    rng = np.random.default_rng(cfg.seed)
    n = geom.n_obs
    G = cfg.n_genes
    xy = np.asarray(geom.obsm["spatial"], dtype=float)
    label, target, is_ntc = _assign_guides(geom, cfg, rng)
    is_pert = (label != "none") & ~is_ntc
    targets = [f"T{k:03d}" for k in range(cfg.n_targets)]
    tidx = {t: k for k, t in enumerate(targets)}

    # truth: autonomous and spillover LFC matrices (targets x genes)
    tau_auto = np.zeros((cfg.n_targets, G))
    m = rng.random((cfg.n_targets, G)) < cfg.frac_auto_nonzero
    tau_auto[m] = rng.normal(0, cfg.auto_lfc_sd, m.sum())
    tau_spill = np.zeros((cfg.n_targets, G))
    spill_targets = rng.random(cfg.n_targets) < cfg.frac_spill_targets
    for k in np.flatnonzero(spill_targets):
        genes = rng.random(G) < cfg.frac_spill_genes
        tau_spill[k, genes] = rng.normal(0, cfg.spill_lfc_sd, genes.sum())

    # baseline expression: gene means (log-normal), cell-type effects, batch, density
    base = rng.lognormal(np.log(cfg.mean_counts_per_cell / G), 1.0, G)
    base = base / base.sum() * cfg.mean_counts_per_cell
    cts = pd.factorize(geom.obs["cell_type"].astype(str))[0]
    ct_eff = rng.normal(0, 0.5, (cts.max() + 1, G))
    samples, s_codes = np.unique(geom.obs["sample"].astype(str).to_numpy(), return_inverse=True)
    batch_eff = (
        rng.normal(0, cfg.batch_sd, (len(samples), G))
        if cfg.batch_sd > 0
        else np.zeros((len(samples), G))
    )
    log_mu = np.log(base)[None, :] + ct_eff[cts] + batch_eff[s_codes]

    # density
    dens = np.zeros(n)
    for s in samples:
        idx = np.flatnonzero(geom.obs["sample"].astype(str).to_numpy() == s)
        tree = cKDTree(xy[idx])
        dens[idx] = np.array([len(x) - 1 for x in tree.query_ball_point(xy[idx], 30.0)])
    log_mu += cfg.beta_density * (dens - dens.mean())[:, None]

    # autonomous effect
    pert_k = np.array([tidx.get(t, -1) for t in target])
    pk = pert_k[is_pert]
    log_mu[is_pert] += tau_auto[pk]

    # spillover: exposure-weighted sum over perturbed neighbours, recipients unperturbed
    spill_dose = np.zeros((n, cfg.n_targets))
    for s in samples:
        idx = np.flatnonzero(geom.obs["sample"].astype(str).to_numpy() == s)
        pidx = idx[is_pert[idx]]
        if len(pidx) == 0:
            continue
        tree = cKDTree(xy[pidx])
        for j in idx:
            nb = tree.query_ball_point(xy[j], cfg.spill_max_um)
            if not nb:
                continue
            nb = np.asarray(nb)
            d = np.linalg.norm(xy[pidx[nb]] - xy[j], axis=1)
            wts = np.exp(-d / cfg.spill_scale_um)
            wts[d == 0] = 0.0
            np.add.at(spill_dose[j], pert_k[pidx[nb]], wts)
    responder = np.ones(n, dtype=bool)
    if cfg.spill_cell_types is not None:
        responder = geom.obs["cell_type"].astype(str).isin(cfg.spill_cell_types).to_numpy()
    spill_lfc = spill_dose @ tau_spill
    spill_lfc[~responder | is_pert] = 0.0
    log_mu += spill_lfc

    mu = np.exp(log_mu)

    # bleed-through on the expected counts: relocate alpha of mu to touching neighbours
    if cfg.alpha_bleed > 0:
        tmp = AnnData(obs=geom.obs.copy())
        tmp.obsm["spatial"] = xy
        build_graph(tmp, GraphConfig(kind="delaunay", max_edge_um=cfg.bleed_max_um))
        A = tmp.obsp["spatial_connectivities"].tocsr()
        deg = np.asarray(A.sum(1)).ravel()
        W = sp.diags(np.where(deg > 0, 1.0 / np.maximum(deg, 1), 0.0)) @ A  # row-normalised
        mu = (1 - cfg.alpha_bleed * (deg > 0)[:, None]) * mu + cfg.alpha_bleed * (W @ mu)

    # negative binomial counts
    phi = cfg.nb_dispersion
    size = 1.0 / phi
    lam = rng.gamma(size, mu / size)
    counts = rng.poisson(lam).astype(np.float32)

    # misassignment: label unassigned neighbours of perturbed cells with the neighbour's guide
    label_obs = label.copy()
    target_obs = target.copy()
    is_ntc_obs = is_ntc.copy()
    n_mis = 0
    if cfg.p_misassign > 0:
        tmp = AnnData(obs=geom.obs.copy())
        tmp.obsm["spatial"] = xy
        build_graph(tmp, GraphConfig(kind="delaunay", max_edge_um=cfg.bleed_max_um))
        A = tmp.obsp["spatial_connectivities"].tocoo()
        src_pert = (label[A.row] != "none") & (label_obs[A.col] == "none")
        cand = np.flatnonzero(src_pert)
        pick = cand[rng.random(len(cand)) < cfg.p_misassign]
        for e in pick:
            j = A.col[e]
            if label_obs[j] == "none":
                label_obs[j] = label[A.row[e]]
                target_obs[j] = target[A.row[e]]
                is_ntc_obs[j] = is_ntc[A.row[e]]
                n_mis += 1
    is_pert_obs = (label_obs != "none") & ~is_ntc_obs

    obs = geom.obs.copy()
    obs["guide"] = label_obs.astype(str)
    obs["guide_confidence"] = np.where(label_obs != "none", 1.0, np.nan)
    obs["target"] = target_obs.astype(str)
    obs["is_ntc"] = is_ntc_obs
    obs["is_perturbed"] = is_pert_obs
    obs["true_guide"] = label.astype(str)
    obs["true_target"] = target.astype(str)
    obs["misassigned"] = label_obs != label
    obs["local_density"] = dens
    adata = AnnData(X=sp.csr_matrix(counts), obs=obs)
    adata.var_names = [f"g{j:03d}" for j in range(G)]
    adata.obsm["spatial"] = xy
    adata.obsm["spill_dose"] = spill_dose
    adata.uns["spatialspill"] = SpatialScreenSchema(
        "simulated", "simulator", "none", notes={"config": str(asdict(cfg))}
    ).to_uns()
    adata.uns["truth"] = {
        "targets": targets,
        "tau_auto": tau_auto,
        "tau_spill": tau_spill,
        "spill_targets": spill_targets,
        "n_misassigned": int(n_mis),
        "config": asdict(cfg),
    }
    return adata


@dataclass
class Scenario:
    """Named simulator settings used in the benchmark grid."""

    name: str
    cfg: SimConfig
    description: str = ""
    tags: list[str] = field(default_factory=list)


def default_scenarios(seed: int = 0) -> list[Scenario]:
    return [
        Scenario(
            "no_effect",
            SimConfig(frac_spill_targets=0.0, frac_auto_nonzero=0.0, seed=seed),
            "no effects at all",
        ),
        Scenario(
            "auto_only",
            SimConfig(frac_spill_targets=0.0, seed=seed),
            "autonomous effects, no spillover",
        ),
        Scenario("spill", SimConfig(seed=seed), "autonomous + spillover"),
        Scenario(
            "spill_bleed",
            SimConfig(alpha_bleed=0.15, seed=seed),
            "+ 15% bleed-through",
            ["artifact"],
        ),
        Scenario(
            "spill_density",
            SimConfig(beta_density=0.05, seed=seed),
            "+ density effect",
            ["artifact"],
        ),
        Scenario(
            "spill_batch", SimConfig(batch_sd=0.3, seed=seed), "+ batch offsets", ["artifact"]
        ),
        Scenario(
            "spill_clonal",
            SimConfig(clonal_radius_um=40.0, clonal_size=6, seed=seed),
            "+ clonal assignment",
            ["design"],
        ),
        Scenario(
            "spill_misassign",
            SimConfig(p_misassign=0.3, seed=seed),
            "+ barcode misassignment",
            ["artifact"],
        ),
        Scenario(
            "everything",
            SimConfig(
                alpha_bleed=0.15,
                beta_density=0.05,
                batch_sd=0.3,
                clonal_radius_um=40.0,
                clonal_size=6,
                p_misassign=0.2,
                seed=seed,
            ),
            "all artifacts",
            ["artifact", "design"],
        ),
    ]
