"""Treated / control group definitions shared by E1 to E3 (docs/estimands.md, Section 3)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from spatialspill.exposure import NTC_PREFIX, Exposure


@dataclass(frozen=True)
class GroupConfig:
    control_policy: str = "ntc"  # "ntc": NTC cells are the reference; "unperturbed": any cell without a targeting guide
    clean_controls: bool = (
        True  # recipients/controls must have no other perturbed neighbour within D_max
    )


def eligible_recipients(
    exp: Exposure, is_ntc: np.ndarray, is_perturbed: np.ndarray, cfg: GroupConfig
) -> np.ndarray:
    if cfg.control_policy == "ntc":
        return is_ntc.copy()
    if cfg.control_policy == "unperturbed":
        return ~is_perturbed
    raise ValueError(cfg.control_policy)


def group_masks(
    exp: Exposure,
    is_ntc: np.ndarray,
    is_perturbed: np.ndarray,
    cfg: GroupConfig,
) -> tuple[sp.csr_matrix, list[sp.csr_matrix], sp.csr_matrix]:
    """Indicator matrices (cells x targets).

    Returns ``(T_auto, [T_ring_b ...], C)`` where ``T_auto[:, k]`` marks cells perturbed with
    target k in an otherwise unexposed neighbourhood, ``T_ring_b[:, k]`` marks eligible recipients
    with target-k neighbours only in ring b, and ``C[:, k]`` marks eligible controls with no
    target-k neighbour (and, if ``clean_controls``, no perturbed neighbour at all).
    """
    n, K = exp.n_cells, len(exp.targets)
    labels = exp.labels
    Z = exp.label_matrix(labels)  # n x K: cell carries target k
    tot = exp.total_count_matrix()  # n x K: number of target-k neighbours within D_max
    any_pert = exp.any_perturbed_within  # n: number of perturbed neighbours (targeting guides)
    elig = eligible_recipients(exp, is_ntc, is_perturbed, cfg)

    # neighbours perturbed with something other than target k
    is_pert_target = np.array(
        [exp.is_perturbed_label.get(t, False) for t in exp.targets], dtype=float
    )
    other_pert = any_pert[:, None] - tot.multiply(is_pert_target[None, :]).toarray()
    clean_other = other_pert <= 0 if cfg.clean_controls else np.ones((n, K), dtype=bool)

    tot_d = tot.toarray()
    T_auto = sp.csr_matrix(Z.multiply(sp.csr_matrix((tot_d == 0) & clean_other)))
    C_d = elig[:, None] & (tot_d == 0) & clean_other
    # NTC pseudo-targets: a cell carrying that very guide is not a control for itself
    for k, t in enumerate(exp.targets):
        if t.startswith(NTC_PREFIX):
            C_d[labels == t, k] = False
    C = sp.csr_matrix(C_d.astype(float))
    rings = []
    for b in range(exp.n_bins):
        cb = exp.counts[b].toarray()
        m = elig[:, None] & (cb >= 1) & (tot_d == cb) & clean_other & (Z.toarray() == 0)
        rings.append(sp.csr_matrix(m.astype(float)))
    return T_auto, rings, C
