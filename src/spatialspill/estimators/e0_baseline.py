"""E0: space-blind and ad hoc baselines used for comparison fairness.

``NeighbourTTest`` reproduces the common published approach to intercellular effects: control
cells with at least one g-perturbed cell within D_max ("with neighbour") are compared with
control cells without one by a Welch t-test per gene, with no strata, no distance bins, no
permutation null and no covariates. The autonomous comparison is perturbed vs control cells by
the same test. This is the "ad hoc statistic" that the ladder E1 to E4 is measured against.

``PseudobulkDE`` is the space-blind autonomous baseline: perturbed and control cells are
aggregated per sample into pseudobulk profiles and compared with a paired t-test over samples
(or a Welch test over cells when only one sample exists). It has no spillover estimand; its
spillover rows are absent, which is the point: a space-blind method cannot report spillover.
"""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy import stats

from spatialspill.estimators.base import EstimateTable, Estimator
from spatialspill.estimators.groups import GroupConfig, eligible_recipients
from spatialspill.exposure import Exposure


def _welch(Ya: np.ndarray, Yb: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ma, mb = Ya.mean(0), Yb.mean(0)
    va, vb = Ya.var(0, ddof=1), Yb.var(0, ddof=1)
    se = np.sqrt(va / len(Ya) + vb / len(Yb))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (ma - mb) / se
        df = (va / len(Ya) + vb / len(Yb)) ** 2 / (
            (va / len(Ya)) ** 2 / (len(Ya) - 1) + (vb / len(Yb)) ** 2 / (len(Yb) - 1)
        )
    p = 2 * stats.t.sf(np.abs(t), np.where(np.isfinite(df), df, 1))
    return ma - mb, se, p


class NeighbourTTest(Estimator):
    name = "E0_neighbour_ttest"

    def __init__(self, min_cells: int = 5, groups: GroupConfig | None = None) -> None:
        self.min_cells = min_cells
        self.groups = groups or GroupConfig()

    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        Y = np.asarray(outcomes, dtype=np.float64)
        obs = adata.obs
        is_ntc = obs["is_ntc"].to_numpy()
        is_pert = obs["is_perturbed"].to_numpy()
        elig = eligible_recipients(exposure, is_ntc, is_pert, self.groups)
        Z = exposure.label_matrix(exposure.labels).tocsc()
        tot = exposure.total_count_matrix().tocsc()
        rows: list[dict[str, object]] = []
        for k, target in enumerate(exposure.targets):
            own = np.asarray(Z[:, k].todense()).ravel() > 0
            tk = np.asarray(tot[:, k].todense()).ravel()
            ctrl = elig & ~own & (tk == 0)
            expo = elig & ~own & (tk >= 1)
            for kind, A, B, ring in (("autonomous", own, ctrl, -1), ("spillover", expo, ctrl, 0)):
                na, nb = int(A.sum()), int(B.sum())
                ident = na >= self.min_cells and nb >= self.min_cells
                if ident:
                    est, se, p = _welch(Y[A], Y[B])
                else:
                    est = se = p = np.full(Y.shape[1], np.nan)
                for g, oname in enumerate(outcome_names):
                    rows.append(
                        {
                            "estimator": self.name,
                            "target": target,
                            "kind": kind,
                            "ring": ring,
                            "cell_type": "all",
                            "outcome": oname,
                            "estimate": float(est[g]),
                            "se": float(se[g]),
                            "ci_low": float(est[g] - 1.96 * se[g]),
                            "ci_high": float(est[g] + 1.96 * se[g]),
                            "pvalue": float(p[g]) if ident and np.isfinite(p[g]) else np.nan,
                            "n_treated": na,
                            "n_control": nb,
                            "identified": bool(ident and np.isfinite(p[g])),
                        }
                    )
        tab = EstimateTable()
        tab.add(rows)
        return tab


class PseudobulkDE(Estimator):
    name = "E0_pseudobulk"

    def __init__(self, min_cells: int = 5, groups: GroupConfig | None = None) -> None:
        self.min_cells = min_cells
        self.groups = groups or GroupConfig()

    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        Y = np.asarray(outcomes, dtype=np.float64)
        obs = adata.obs
        is_ntc = obs["is_ntc"].to_numpy()
        is_pert = obs["is_perturbed"].to_numpy()
        elig = eligible_recipients(exposure, is_ntc, is_pert, self.groups)
        samples = obs["sample"].astype(str).to_numpy()
        uniq = np.unique(samples)
        Z = exposure.label_matrix(exposure.labels).tocsc()
        rows: list[dict[str, object]] = []
        for k, target in enumerate(exposure.targets):
            own = np.asarray(Z[:, k].todense()).ravel() > 0
            ctrl = elig & ~own
            na, nb = int(own.sum()), int(ctrl.sum())
            ident = na >= self.min_cells and nb >= self.min_cells
            if ident and len(uniq) >= 3:
                diffs = []
                for s in uniq:
                    m = samples == s
                    if (own & m).sum() >= 2 and (ctrl & m).sum() >= 2:
                        diffs.append(Y[own & m].mean(0) - Y[ctrl & m].mean(0))
                if len(diffs) >= 3:
                    D = np.stack(diffs)
                    est = D.mean(0)
                    se = D.std(0, ddof=1) / np.sqrt(len(D))
                    with np.errstate(divide="ignore", invalid="ignore"):
                        p = 2 * stats.t.sf(np.abs(est / se), len(D) - 1)
                else:
                    est, se, p = _welch(Y[own], Y[ctrl])
            elif ident:
                est, se, p = _welch(Y[own], Y[ctrl])
            else:
                est = se = p = np.full(Y.shape[1], np.nan)
            for g, oname in enumerate(outcome_names):
                rows.append(
                    {
                        "estimator": self.name,
                        "target": target,
                        "kind": "autonomous",
                        "ring": -1,
                        "cell_type": "all",
                        "outcome": oname,
                        "estimate": float(est[g]),
                        "se": float(se[g]),
                        "ci_low": float(est[g] - 1.96 * se[g]),
                        "ci_high": float(est[g] + 1.96 * se[g]),
                        "pvalue": float(p[g]) if ident and np.isfinite(p[g]) else np.nan,
                        "n_treated": na,
                        "n_control": nb,
                        "identified": bool(ident and np.isfinite(p[g])),
                    }
                )
        tab = EstimateTable()
        tab.add(rows)
        return tab
