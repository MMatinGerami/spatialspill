"""Dataset loaders producing the spatialspill AnnData schema (see spatialspill.schema)."""

from spatialspill.loaders.perturb_dbit import load_perturb_dbit
from spatialspill.loaders.perturb_map import load_perturb_map
from spatialspill.loaders.perturb_multi import load_perturb_multi

__all__ = ["load_perturb_dbit", "load_perturb_map", "load_perturb_multi"]
