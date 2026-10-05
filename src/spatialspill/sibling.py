"""Undetected-sibling model: recover spillover from unassigned recipients.

Under clonal growth or local delivery with incomplete guide calling, an unassigned cell next to
a g cell is, with probability pi_g, an undetected g cell (a sibling) and otherwise a true
recipient. Its expected profile (relative to unassigned cells away from any g cell) is then

    u_g = pi_g * a_g + (1 - pi_g) * s_g

where a_g is g's autonomous profile and s_g the spillover profile onto true recipients. With
a_g estimated from detected g cells and u_g from unassigned neighbours, pi_g is estimated by
least squares projection of u_g on a_g across outcomes (clipped to [0, 1]) and the corrected
spillover is s_hat = (u_g - pi_hat a_g) / (1 - pi_hat). Uncertainty: bootstrap over the
unassigned neighbours (percentile intervals). Identification needs ||a_g|| well above noise and
assumes true recipients do not reproduce g's autonomous profile (the same assumption the
bleed-through detector makes; a biological "same programme" spillover is absorbed into pi).

This is a mixture deconvolution, not a cell-level classifier. It is tested on the simulator's
clonal and misassignment scenarios (scripts/benchmark_sibling.py) before any real-data use.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from anndata import AnnData

from spatialspill.exposure import Exposure


@dataclass
class SiblingEstimate:
    target: str
    pi: float
    pi_ci: tuple[float, float]
    n_unassigned_neighbours: int
    n_target_cells: int
    corrected: np.ndarray  # spillover profile over outcomes
    corrected_ci: np.ndarray  # (2, n_outcomes)
    uncorrected: np.ndarray
    autonomous: np.ndarray


def _profile(Y: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return Y[mask].mean(axis=0)


def estimate_pi(u: np.ndarray, a: np.ndarray) -> float:
    denom = float(a @ a)
    if denom <= 0:
        return float("nan")
    return float(np.clip((u @ a) / denom, 0.0, 1.0))


def sibling_corrected_spillover(
    adata: AnnData,
    exposure: Exposure,
    Y: np.ndarray,
    target: str,
    ring: int | None = None,
    reference: str = "ntc",
    n_boot: int = 200,
    min_cells: int = 10,
    seed: int = 0,
) -> SiblingEstimate | None:
    """Mixture-corrected spillover of ``target`` onto unassigned neighbours.

    ``ring`` restricts exposure to one ring (None: any ring within D_max). ``reference`` is the
    group that defines zero: "ntc" (NTC cells with no g neighbour) or "unassigned" (unassigned
    cells with no perturbed neighbour at all).
    """
    k = exposure.targets.index(target)
    obs = adata.obs
    labels = exposure.labels
    own = labels == target
    unassigned = np.asarray(labels == "none")
    counts = exposure.counts[ring] if ring is not None else exposure.total_count_matrix()
    c = np.asarray(counts[:, k].todense()).ravel()
    tot = np.asarray(exposure.total_count_matrix()[:, k].todense()).ravel()
    nb = unassigned & (c > 0)
    if reference == "ntc":
        ref = obs["is_ntc"].to_numpy() & (tot == 0)
    else:
        ref = unassigned & (exposure.any_perturbed_within == 0)
    if nb.sum() < min_cells or own.sum() < min_cells or ref.sum() < min_cells:
        return None
    base = _profile(Y, ref)
    a = _profile(Y, own) - base
    u = _profile(Y, nb) - base
    pi = estimate_pi(u, a)
    corr = (u - pi * a) / max(1 - pi, 1e-3)
    rng = np.random.default_rng(seed)
    nb_idx = np.flatnonzero(nb)
    own_idx = np.flatnonzero(own)
    boots = np.empty((n_boot, Y.shape[1]))
    pis = np.empty(n_boot)
    for b in range(n_boot):
        ub = Y[rng.choice(nb_idx, len(nb_idx))].mean(0) - base
        ab = Y[rng.choice(own_idx, len(own_idx))].mean(0) - base
        pb = estimate_pi(ub, ab)
        pis[b] = pb
        boots[b] = (ub - pb * ab) / max(1 - pb, 1e-3)
    return SiblingEstimate(
        target=target,
        pi=pi,
        pi_ci=(float(np.nanpercentile(pis, 2.5)), float(np.nanpercentile(pis, 97.5))),
        n_unassigned_neighbours=int(nb.sum()),
        n_target_cells=int(own.sum()),
        corrected=corr,
        corrected_ci=np.nanpercentile(boots, [2.5, 97.5], axis=0),
        uncorrected=u,
        autonomous=a,
    )


def sibling_table(
    adata: AnnData,
    exposure: Exposure,
    Y: np.ndarray,
    outcome_names: list[str],
    ring: int | None = None,
    **kw: object,
) -> pd.DataFrame:
    """One row per (target, outcome): uncorrected and corrected spillover with intervals, plus pi."""
    rows = []
    for t in exposure.targets:
        if t.startswith("NTC:"):
            continue
        est = sibling_corrected_spillover(adata, exposure, Y, t, ring=ring, **kw)  # type: ignore[arg-type]
        if est is None:
            continue
        for g, name in enumerate(outcome_names):
            rows.append(
                {
                    "target": t,
                    "outcome": name,
                    "ring": -1 if ring is None else ring,
                    "pi": est.pi,
                    "pi_low": est.pi_ci[0],
                    "pi_high": est.pi_ci[1],
                    "n_unassigned_neighbours": est.n_unassigned_neighbours,
                    "uncorrected": float(est.uncorrected[g]),
                    "corrected": float(est.corrected[g]),
                    "corrected_low": float(est.corrected_ci[0, g]),
                    "corrected_high": float(est.corrected_ci[1, g]),
                    "autonomous": float(est.autonomous[g]),
                }
            )
    return pd.DataFrame(rows)
