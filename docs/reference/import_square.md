# Square Analysis Import

## Overview

`import_visiumhd_square` packages Space Ranger square-bin analysis outputs (clusters, DE, and optional UMAP coordinates) into PMTiles and a square assets JSON for CartLoader. Use it directly, or via `run_visiumhd --import-squares` to fold square-bin layers into your catalog.

---
## Action

Convert one bin size of Visium HD square results into PMTiles, write `<square-id>-sqXXX_assets.json`, and optionally update an existing `catalog.yaml`.

---
## Requirements

- Space Ranger square-bin outputs for a given bin size (8 µm or 16 µm are common):
  - Positions Parquet (`spatial/tissue_positions.parquet`)
  - Cluster CSV (`analysis/clustering/.../clusters.csv`)
  - Differential expression CSV (`analysis/diffexp/.../differential_expression.csv`)
  - (Optional) UMAP CSV (`analysis/pca/gene_expression_10_components/projection.csv`)
- Scale factors JSON (`spatial/scalefactors_json.json`) to derive units-per-µm, or provide `--units-per-um`.
- Tippecanoe and parquet-tools (or pigz/polars) available on PATH.

---
## Example Usage

=== "JSON mode (preferred; auto paths)"
    ```bash
    cartloader import_visiumhd_square \
      --in-json /path/to/space_ranger_assets.json \
      --bin-size 8 \
      --outprefix /path/to/cartload2/squares/spaceranger-sq008 \
      --threads 8
    ```

=== "Manual mode (explicit files)"
    ```bash
    cartloader import_visiumhd_square \
      --in-dir /path/to/space_ranger_output \
      --bin-size 8 \
      --parquet spatial/tissue_positions.parquet \
      --csv-clust analysis/clustering/gene_expression_graphclust/clusters.csv \
      --csv-diffexp analysis/diffexp/gene_expression_graphclust/differential_expression.csv \
      --csv-umap analysis/pca/gene_expression_10_components/projection.csv \
      --outprefix /path/to/cartload2/squares/spaceranger-sq008 \
      --threads 8
    ```

---
## Parameters

### Input/Output
- `--in-json` (str): Space Ranger assets JSON (from `load_space_ranger` or `run_visiumhd --load-space-ranger`); mutually exclusive with manual mode.
- `--in-dir` (str): Space Ranger output directory for manual mode.
- `--bin-size` (int, required): Bin size in µm (e.g., 8 or 16).
- `--outprefix` (str, required): Output prefix for generated artifacts.
- `--id`, `--name` (str): Override asset ID/display name (default: basename of `--outprefix`).

### Manual Inputs (relative to `--in-dir`; override as needed)
- `--parquet` (str): Positions Parquet (default: `spatial/tissue_positions.parquet`).
- `--csv-clust` (str): Cluster assignments (default: `analysis/clustering/gene_expression_graphclust/clusters.csv`).
- `--csv-diffexp` (str): Differential expression (default: `analysis/diffexp/gene_expression_graphclust/differential_expression.csv`).
- `--csv-umap` (str): UMAP projection CSV (default: `analysis/pca/gene_expression_10_components/projection.csv`; optional).
- `--scale-json` (str): Scale factors JSON (default: `spatial/scalefactors_json.json`).
- `--units-per-um` (float): Override microns-per-unit if scale JSON is missing.

### Column Overrides
- `--pos-colname-*`, `--clust-colname-*`, `--umap-colname-*`: Adjust barcode/coordinate column names if they differ from defaults.

### Conversion / Catalog
- `--min-zoom`, `--max-zoom`: Tile zooms for square PMTiles.
- `--umap-min-zoom`, `--umap-max-zoom`: Tile zooms for UMAP PMTiles.
- `--max-tile-bytes`, `--max-feature-counts`, `--preserve-point-density-thres`: Tippecanoe tuning knobs.
- `--update-catalog` (flag): Append square assets to an existing `catalog.yaml` (default path `<in-dir>/catalog.yaml` or `--catalog-yaml`).
- `--catalog-yaml` (str): Catalog path to update when `--update-catalog` is set.

### Run / Env
- `--threads` (int): Threads for tippecanoe.
- `--use-parquet-tools` (flag): Use `parquet-tools csv` instead of polars/pigz conversion.
- `--tsv-cmap` (str): Cluster color map TSV (default: `assets/fixed_color_map_60.tsv`).
- `--tippecanoe`, `--parquet-tools`, `--pigz`, `--gzip`, `--R`: Override tool paths.
- `--log`, `--log-suffix`: Write a log next to outputs.
- `--restart`, `--dry-run`: Inherit from wrappers (when called via `run_visiumhd`).

---
## Outputs

Written next to `--outprefix`:

- `<outprefix>.pmtiles`: Square-bin clusters as tiles (one layer per bin size).
- `<outprefix>-bulk-de.tsv`: Differential expression per cluster.
- `<outprefix>-rgb.tsv`: Cluster colormap TSV.
- `<outprefix>-umap.tsv.gz` / `<outprefix>-umap.pmtiles` / `<outprefix>.umap.png` / `<outprefix>.umap.single.binary.png`: UMAP coordinates and visualizations (if `--csv-umap` provided).
- `<outprefix>_assets.json`: Square assets manifest referenced by `catalog.yaml`.

If `--update-catalog` is set, the manifest is appended to the existing catalog under the provided `catalog.yaml`.
