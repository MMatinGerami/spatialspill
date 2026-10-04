"""Held-out-gene prediction of perturbation effect profiles from gene embeddings (C5).

Given an estimate table (E1/E2 output) and one or more embedding matrices, build for each
target a response vector (the z-statistics of its autonomous or spillover effects over
outcomes), then predict the response of a held-out target from the responses of the other
targets through their embeddings (kernel ridge on the embedding space; leave-one-target-out).
Scores: Pearson r between predicted and observed response, per target and overall, against a
null that permutes the embedding rows (same model, meaningless features). A failure to beat
the null is reported as a finding.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.kernel_ridge import KernelRidge


def response_matrix(
    estimates: pd.DataFrame,
    kind: str,
    cell_type: str = "all",
    rings: list[int] | None = None,
    stat: str = "z",
) -> pd.DataFrame:
    """Targets x outcomes matrix of z = estimate / se (or estimates) for one effect kind.

    Spillover rings are averaged (inverse-variance weighted) when several are given.
    """
    sub = estimates[
        (estimates["kind"] == kind)
        & (estimates["cell_type"] == cell_type)
        & estimates["identified"].astype(bool)
    ]
    if kind == "spillover" and rings is not None:
        sub = sub[sub["ring"].isin(rings)]
    sub = sub[~sub["target"].astype(str).str.startswith("NTC:")]
    if sub.empty:
        return pd.DataFrame()
    sub = sub.assign(w=1.0 / sub["se"].replace(0, np.nan) ** 2)
    if stat == "z":
        sub = sub.assign(v=sub["estimate"] / sub["se"].replace(0, np.nan))
    else:
        sub = sub.assign(v=sub["estimate"])
    sub = sub.assign(vw=sub["v"] * sub["w"])
    g = sub.groupby(["target", "outcome"])[["vw", "w"]].sum(min_count=1)
    num = g["vw"] / g["w"].replace(0, np.nan)
    M = num.unstack("outcome")
    return M.dropna(axis=0, how="all")


def loto_predict(
    R: pd.DataFrame,
    E: pd.DataFrame,
    alpha: float = 1.0,
    gamma: float | None = None,
    seed: int = 0,
    n_null: int = 20,
) -> dict[str, object]:
    """Leave-one-target-out kernel ridge prediction of response rows from embeddings.

    Returns per-target Pearson r, the mean r, and the mean r under row-permuted embeddings.
    """
    common = R.index.intersection(E.index)
    if len(common) < 5:
        return {
            "n_targets": len(common),
            "mean_r": np.nan,
            "null_mean_r": np.nan,
            "per_target": pd.DataFrame(),
        }
    Rm = R.loc[common].to_numpy(dtype=float)
    Rm = np.where(np.isfinite(Rm), Rm, 0.0)
    Em = E.loc[common].to_numpy(dtype=float)
    Em = (Em - Em.mean(0)) / (Em.std(0) + 1e-9)
    if gamma is None:
        gamma = 1.0 / Em.shape[1]

    def run(Emat: np.ndarray) -> np.ndarray:
        rs = np.empty(len(common))
        for i in range(len(common)):
            tr = np.arange(len(common)) != i
            model = KernelRidge(alpha=alpha, kernel="rbf", gamma=gamma).fit(Emat[tr], Rm[tr])
            pred = model.predict(Emat[i : i + 1])[0]
            obs = Rm[i]
            rs[i] = np.corrcoef(pred, obs)[0, 1] if np.std(pred) > 0 and np.std(obs) > 0 else 0.0
        return rs

    rs = run(Em)
    rng = np.random.default_rng(seed)
    null = np.array([run(Em[rng.permutation(len(common))]).mean() for _ in range(n_null)])
    per = pd.DataFrame({"target": common, "r": rs})
    return {
        "n_targets": len(common),
        "mean_r": float(np.mean(rs)),
        "null_mean_r": float(null.mean()),
        "null_sd_r": float(null.std()),
        "z_vs_null": float((np.mean(rs) - null.mean()) / null.std()) if null.std() > 0 else np.nan,
        "per_target": per,
    }
