# Square Analysis Import

## Overview

`import_visiumhd_square` packages Space Ranger square-bin analysis outputs (clusters, DE, and optional UMAP coordinates) into PMTiles and a square assets JSON for CartLoader. Use it directly, or via `run_visiumhd --import-squares` to fold square-bin layers into your catalog.

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
## Actions

Convert one bin size of Visium HD square results into PMTiles, write `<square-id>-sqXXX_assets.json`, and optionally update an existing `catalog.yaml`.

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
- `--outprefix` (str, required): Output prefix for generated artifacts.
- `--id`, `--name` (str): Override asset ID/display name (default: basename of `--outprefix`).

### Manual Inputs (relative to `--in-dir`; override as needed)
- `--in-dir` (str): Space Ranger output directory for manual mode.
- `--bin-size` (int, required): Bin size in µm (e.g., 8 or 16).
- `--parquet` (str): Positions Parquet (default: `spatial/tissue_positions.parquet`).
- `--csv-clust` (str): Cluster assignments (default: `analysis/clustering/gene_expression_graphclust/clusters.csv`).
- `--csv-diffexp` (str): Differential expression (default: `analysis/diffexp/gene_expression_graphclust/differential_expression.csv`).
- `--csv-umap` (str): UMAP projection CSV (default: `analysis/pca/gene_expression_10_components/projection.csv`; optional).
- `--scale-json` (str): Scale factors JSON (default: `spatial/scalefactors_json.json`).
- `--units-per-um` (float): Override microns-per-unit if scale JSON is missing.

### Column Overrides
- `--pos-colname-*`, `--clust-colname-*`, `--umap-colname-*`: Adjust barcode/coordinate column names if they differ from defaults.


??? note "Auxiliary Parameters"
    Recommend to use the default values; override only if needed.

    **PMTiles Conversion**:
    * `--min-zoom` (int): Minimum zoom level (default: 0).
    * `--max-zoom` (int): Maximum zoom level (default: 18).
    * `--max-tile-bytes` (int): Max bytes per tile.
    * `--preserve-point-density-thres` (int): Density preservation threshold.
    * `--max-feature-counts` (int): Max features per tile.

    **Auxiliary UMAP Parameters**:
    * `--umap-colname-factor` (str): Column name for the clustering factor (default: `Cluster`).
    * `--umap-colname-x` (str): Column name for UMAP X (default: `UMAP1`).
    * `--umap-colname-y` (str): Column name for UMAP Y (default: `UMAP2`).
    * `--umap-min-zoom` (int): Min UMAP zoom (default: 0).
    * `--umap-max-zoom` (int): Max UMAP zoom (default: 18).

    **Auxiliary Colname Parameters**:
    * `--pos-colname-barcode` (str): Column name for barcode in -- (default: `Barcode`).
    * `--pos-colname-x` (str): Column name for X position.
    * `--pos-colname-y` (str): Column name for Y position.
    * `--cluster-colname-barcode` (str): Column name for the barcode (default: `Barcode`).
    * `--cluster-colname-cluster` (str): Column name for cluster ID (default: `Cluster`).

    **Auxiliary Catalog Parameters**:
    * `--update-catalog` (flag): Append square assets to an existing `catalog.yaml` (default path `<in-dir>/catalog.yaml` or `--catalog-yaml`).
    * `--catalog-yaml` (str): Catalog path to update when `--update-catalog` is set.

    **Auxiliary Run Parameters**:
    * `--threads` (int): Threads for tippecanoe.
    - `--tsv-cmap` (str): Cluster color map TSV (default: `assets/fixed_color_map_60.tsv`).
    - `--log`, `--log-suffix`: Write a log next to outputs.

    **Auxiliary Environment Parameters**:
    - `--use-parquet-tools` (flag): Use `parquet-tools csv` instead of polars/pigz conversion.
    * `--parquet-tools` (str): Path to `parquet-tools` binary (default: `parquet-tools`).
    * `--tippecanoe` (str): Path to `tippecanoe` binary (default: `tippecanoe`).
    * `--pigz` (str): Path to `pigz` binary (default: `pigz`).
    * `--pigz-threads` (int): Threads for pigz (default: 4).
    * `--gzip` (str): Path to `gzip` (default: `gzip`).

---
## Output

Written next to `--outprefix`:

- `<outprefix>.pmtiles`: Square-bin clusters as tiles (one layer per bin size).
- `<outprefix>-bulk-de.tsv`: Differential expression per cluster.
- `<outprefix>-rgb.tsv`: Cluster colormap TSV.
- `<outprefix>-umap.tsv.gz` / `<outprefix>-umap.pmtiles` / `<outprefix>.umap.png` / `<outprefix>.umap.single.binary.png`: UMAP coordinates and visualizations (if `--csv-umap` provided).
- `<outprefix>_assets.json`: Square assets manifest referenced by `catalog.yaml`.

If `--update-catalog` is set, the manifest is appended to the existing catalog under the provided `catalog.yaml`.
