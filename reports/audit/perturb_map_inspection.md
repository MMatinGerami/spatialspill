# Perturb-map (GSE193460) raw data inspection

Dhainaut et al., Cell 2022. Four Visium sections of KP lung tumours with Pro-Code CRISPR
knockouts, Space Ranger outputs in `data/raw/perturb_map/extracted/`. Numbers below were
produced by `load_perturb_map("data/raw/perturb_map")` plus the snippet at the end; the loader
is `src/spatialspill/loaders/perturb_map.py`. `validate(adata)` returns no problems; the merged
object is 6363 spots x 32289 genes.

## Per section

| section | mouse | spots kept | off_tissue dropped | tumor (kmeans) | labelled lesion spots | periphery spots | knockouts labelled | median counts | median genes | px per um (pitch) |
|:--------|------:|-----------:|-------------------:|---------------:|----------------------:|----------------:|:-------------------|--------------:|-------------:|------------------:|
| KP_1    | 4.1   | 1903       | 3                  | 540            | 329                   | 149             | Jak2, Tgfbr2       | 8005          | 3284         | 0.4901            |
| KP_2    | 4.3   | 1872       | 59                 | 639            | 383                   | 194             | Ifngr2, Tgfbr2     | 8805          | 3235         | 0.4901            |
| KP_3    | 5.1   | 1233       | 11                 | 713            | 350                   | 148             | Ifngr2             | 10469         | 3632         | 0.4901            |
| KP_4    | 7.2   | 1355       | 1191               | 724            | 349                   | 104             | Tgfbr2             | 10375         | 3460         | 0.4901            |

"spots kept" = `in_tissue == 1` spots in the filtered matrix minus the authors' `off_tissue`
list. After this filter the barcodes equal the rows of `spot_annotation.csv.gz` exactly in all
four sections. "labelled lesion spots" are spots with a phenotype other than NA or `periphery`.

## Spots per phenotype label

| sample | guide (phenotype) | spots | target       |
|:-------|:------------------|------:|:-------------|
| KP_1   | Jak2_1            | 21    | Jak2         |
| KP_1   | KP_1-1            | 110   | KP_unlabeled |
| KP_1   | KP_1-2            | 75    | KP_unlabeled |
| KP_1   | KP_1-3            | 48    | KP_unlabeled |
| KP_1   | Tgfbr2_1          | 75    | Tgfbr2       |
| KP_1   | none              | 1425  | none         |
| KP_1   | periphery         | 149   | none         |
| KP_2   | Ifngr2_2          | 41    | Ifngr2       |
| KP_2   | KP_2-1            | 148   | KP_unlabeled |
| KP_2   | KP_2-2            | 97    | KP_unlabeled |
| KP_2   | KP_2-3            | 35    | KP_unlabeled |
| KP_2   | Tgfbr2_2          | 62    | Tgfbr2       |
| KP_2   | none              | 1295  | none         |
| KP_2   | periphery         | 194   | none         |
| KP_3   | Ifngr2_3          | 64    | Ifngr2       |
| KP_3   | KP_3-1            | 133   | KP_unlabeled |
| KP_3   | KP_3-2            | 98    | KP_unlabeled |
| KP_3   | KP_3-3            | 55    | KP_unlabeled |
| KP_3   | none              | 735   | none         |
| KP_3   | periphery         | 148   | none         |
| KP_4   | KP_4-1            | 115   | KP_unlabeled |
| KP_4   | KP_4-2            | 74    | KP_unlabeled |
| KP_4   | KP_4-3            | 38    | KP_unlabeled |
| KP_4   | KP_4-4            | 35    | KP_unlabeled |
| KP_4   | Tgfbr2_4-1        | 50    | Tgfbr2       |
| KP_4   | Tgfbr2_4-2        | 37    | Tgfbr2       |
| KP_4   | none              | 902   | none         |
| KP_4   | periphery         | 104   | none         |

Named knockouts across the dataset: Tgfbr2 (224 spots in 4 sections), Ifngr2 (105 spots in
2 sections), Jak2 (21 spots in 1 section). Every labelled lesion cluster is one Leiden cluster
(`leiden_clusters` 1 to 7 per section, 1 is always `periphery`). No non-targeting control label
exists, so `is_ntc` is False everywhere.

`cell_type` combines the authors' k-means call with the lesion status:

| cell_type        | spots |
|:-----------------|------:|
| normal_unlabeled | 3548  |
| tumor_lesion     | 1389  |
| tumor_unlabeled  | 809   |
| tumor_periphery  | 418   |
| normal_periphery | 177   |
| normal_lesion    | 22    |

## Notes and surprises

- `tissue_positions_list.csv.gz` has no header (Space Ranger v1 layout). Columns are barcode,
  in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres.
- The "fullres" pixel coordinates are on a downsampled image: the median nearest-neighbour
  distance between spot centres is 49.0 px in every section, so 100 um = 49.0 px
  (0.490 px/um, pitch estimated over all 4992 array positions).
- `spot_diameter_fullres` is 32.32 px in every section. Taken as 55 um it implies 0.588 px/um,
  which would put the array pitch at 84 um instead of 100 um. The ratio diameter / pitch is
  0.66, so the Space Ranger diameter corresponds to about 65 um (a known property of the
  scalefactor). The loader therefore anchors on the 100 um pitch by default and exposes
  `micron_scale="spot_diameter"` for comparison; both scales are stored in
  `uns["spatialspill"]["notes"]`.
- KP_4 has 1191 off-tissue spots (47 percent of its filtered matrix). These are low-count spots
  (median counts over the filtered matrix 3896 versus 10375 after removal), so the off-tissue
  list matters for this section.
- 40 gene symbols are duplicated in the Space Ranger feature list (32289 features, reference
  `refdata-gex-mm10-2020-A-mCherry`, which includes an mCherry transgene);
  `var_names_make_unique()` is applied and Ensembl ids are kept in `var["gene_ids"]`.
- `nCount_Spatial` in the annotation equals the row sum of the filtered matrix (same median, to
  within rounding of the medians across different spot sets), so X is unnormalised.
- Edge distance uses the bounding-box proxy from `schema.edge_distance_from_coords`; no tissue
  mask is parsed from the images.

## Snippet used

Run from the repo root with `uv run --with tabulate python <file>` (tabulate is only needed for
`to_markdown`).

```python
import pandas as pd
from spatialspill.loaders.perturb_map import load_perturb_map
from spatialspill.schema import validate

adata = load_perturb_map("data/raw/perturb_map")
print("validate:", validate(adata), adata.shape)
rows = []
for s, o in adata.obs.groupby("sample", observed=True):
    named = sorted(o.loc[o["is_perturbed"], "target"].unique())
    rows.append(
        {
            "section": s,
            "mouse": o["batch"].iloc[0].removeprefix("mouse_"),
            "spots_kept": len(o),
            "off_tissue_dropped": adata.uns["spatialspill"]["notes"][f"{s}_off_tissue_dropped"],
            "tumor_kmeans": int(o["cell_type"].str.startswith("tumor").sum()),
            "labelled_lesion_spots": int((~o["guide"].isin(["none", "periphery"])).sum()),
            "periphery_spots": int((o["guide"] == "periphery").sum()),
            "knockouts": ", ".join(named),
            "median_counts": int(o["n_counts"].median()),
            "median_genes": int(o["n_features"].median()),
            "px_per_um_pitch": adata.uns["spatialspill"]["notes"][f"{s}_px_per_um_pitch"],
        }
    )
print(pd.DataFrame(rows).to_markdown(index=False))
phen = adata.obs.groupby(["sample", "guide"], observed=True).size().reset_index(name="spots")
phen["target"] = (
    adata.obs.drop_duplicates("guide").set_index("guide").loc[phen["guide"], "target"].to_numpy()
)
print(phen.to_markdown(index=False))
print(adata.obs["cell_type"].value_counts().to_markdown())
```

The raw-file checks (header detection, barcode set comparisons, scalefactor values, duplicate
gene names) were done with `pandas.read_csv`, `gzip`/`json` and `scanpy.read_10x_h5` on each
section before the loader was written; the loader encodes the same logic.
