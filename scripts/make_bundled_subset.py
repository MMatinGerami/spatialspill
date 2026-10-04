"""Create the tiny bundled subset used by CI and `make smoke`.

Perturb-DBiT sample mTSG.LV.2 (2,500 pixels), restricted to the 300 most-detected genes.
Deterministic; re-running produces an identical file.
"""

from __future__ import annotations

from pathlib import Path

from spatialspill.loaders.perturb_dbit import load_perturb_dbit
from spatialspill.schema import validate

OUT = (
    Path(__file__).resolve().parents[1] / "data" / "bundled" / "perturb_dbit_mTSG.LV.2_top300.h5ad"
)


def main() -> None:
    a = load_perturb_dbit(samples=["mTSG.LV.2"], max_genes=300)
    validate(a)
    a.uns["spatialspill"]["notes"]["bundled"] = "top 300 genes by total counts; for CI only"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(OUT, compression="gzip")
    print(a, OUT.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
