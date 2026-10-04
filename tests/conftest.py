from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from spatialspill.schema import SpatialScreenSchema, edge_distance_from_coords


def make_toy(
    n: int = 300,
    n_genes: int = 20,
    n_samples: int = 2,
    seed: int = 0,
    field_um: float = 600.0,
    p_perturbed: float = 0.2,
) -> AnnData:
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, field_um, size=(n, 2))
    sample = np.array([f"s{i % n_samples}" for i in range(n)])
    pp = p_perturbed / 2
    targets = rng.choice(
        ["NTC", "GENE_A", "GENE_B", "none"], size=n, p=[0.7 - p_perturbed, pp, pp, 0.3]
    )
    obs = pd.DataFrame(
        {
            "guide": [f"{t}_g{rng.integers(1, 4)}" if t != "none" else "none" for t in targets],
            "guide_confidence": np.where(targets == "none", np.nan, rng.uniform(0.5, 1, n)),
            "target": targets,
            "is_ntc": targets == "NTC",
            "is_perturbed": ~np.isin(targets, ["NTC", "none"]),
            "area": rng.uniform(50, 150, n),
            "cell_type": rng.choice(["A", "B"], size=n),
            "sample": sample,
            "batch": sample,
            "edge_distance": edge_distance_from_coords(xy, sample),
        },
        index=[f"c{i}" for i in range(n)],
    )
    X = rng.poisson(2.0, size=(n, n_genes)).astype(np.float32)
    adata = AnnData(X=X, obs=obs)
    adata.var_names = [f"G{j}" for j in range(n_genes)]
    adata.obsm["spatial"] = xy
    adata.uns["spatialspill"] = SpatialScreenSchema("toy", "toy", "none").to_uns()
    return adata


@pytest.fixture
def toy() -> AnnData:
    return make_toy()
