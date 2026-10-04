"""Shared result container and interface for estimators."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from anndata import AnnData

from spatialspill.exposure import Exposure

COLUMNS = [
    "estimator",
    "target",
    "kind",  # "autonomous" | "spillover"
    "ring",  # -1 for autonomous
    "cell_type",
    "outcome",
    "estimate",
    "se",
    "ci_low",
    "ci_high",
    "pvalue",
    "n_treated",
    "n_control",
    "identified",
]


@dataclass
class EstimateTable:
    """Long table of estimates; one row per (target, kind, ring, cell type, outcome)."""

    df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=COLUMNS))

    def add(self, rows: list[dict[str, object]]) -> None:
        if rows:
            new = pd.DataFrame(rows)
            self.df = pd.concat([self.df, new], ignore_index=True) if len(self.df) else new

    def with_fdr(self, q_col: str = "qvalue") -> pd.DataFrame:
        from spatialspill.multiple_testing import bh_fdr

        df = self.df.copy()
        df[q_col] = np.nan
        m = df["identified"].astype(bool).to_numpy()
        df.loc[m, q_col] = bh_fdr(df.loc[m, "pvalue"].to_numpy(dtype=float))
        return df


class Estimator(ABC):
    name: str = "base"

    @abstractmethod
    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        """Estimate autonomous and spillover effects for every target in ``exposure``.

        ``outcomes`` is (n_cells x n_outcomes) dense array on the analysis scale.
        """
