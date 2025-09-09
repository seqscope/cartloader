# Cell Analysis Import

## Overview

Import platform‑specific cell analysis outputs (cells, boundaries, clusters/DE), convert them to web‑ready layers (PMTiles/GeoJSON), and optionally register them in a cartloader catalog. Supported platforms include 10x Xenium and 10x Visium HD.

## Requirements

- Cell analysis outputs from the platform (see per‑platform sections).
- For PMTiles conversion: `tippecanoe` installed and on PATH.
- Optional catalog update: an existing `catalog.yaml` from [Asset Packaging](./run_cartload2.md).

## Actions

- `--cells`: Convert cell points/centroids to PMTiles with counts and top cluster index (topK).
- `--boundaries`: Convert cell polygons to GeoJSON with `cell_id` and `topK` attributes.
- `--summary`: Write a JSON summary of inputs/outputs and basic stats.
- `--update-catalog`: Append generated layers to an existing `catalog.yaml` (not included in `--all`).
- `--all`: Enable `--cells`, `--boundaries`, and `--summary` together.

## Example Usages

CartLoader accepts two input modes for each platform:

- JSON mode: pass `--in-json` with required keys.
- Manual mode: pass `--indir` with specific file locations.

=== "Xenium Example"

    ```bash
    # JSON input mode (requires keys: CELL, BOUNDARY, CLUSTER, DE)
    cartloader import_xenium_cell \
      --cells --boundaries --summary \
      --outprefix /path/to/out/xeniumranger \
      --in-json /path/to/xenium_cells.json \
      --tippecanoe /path/to/tippecanoe

    # Manual input mode
    cartloader import_xenium_cell \
      --cells --boundaries --summary \
      --outprefix /path/to/out/xeniumranger \
      --indir /path/to/xenium_ranger_output \
      --csv-cells analysis/cell_feature_matrix/cells.csv \
      --csv-boundaries analysis/cell_segmentation/cells_boundary.geojson \
      --csv-clust analysis/clustering/gene_expression_graphclust/clusters.csv \
      --csv-diffexp analysis/diffexp/gene_expression_graphclust/differential_expression.csv \
      --tippecanoe /path/to/tippecanoe
    ```

=== "Visium HD Example"

    ```bash
    # JSON input mode (requires keys: CELL_FEATURE_MEX, CELL_GEOJSON, CLUSTER, DE)
    cartloader import_visiumhd_cell \
      --cells --boundaries --summary \
      --outprefix /path/to/out/spaceranger \
      --in-json /path/to/visiumhd_cells.json

    # Manual input mode
    cartloader import_visiumhd_cell \
      --cells --boundaries --summary \
      --outprefix /path/to/out/spaceranger \
      --indir /path/to/spaceranger_output \
      --mtx-cell filtered_feature_bc_matrix \
      --geojson-cells spatial/cells.geojson \
      --csv-clust analysis/clustering/gene_expression_graphclust/clusters.csv \
      --csv-diffexp analysis/diffexp/gene_expression_graphclust/differential_expression.csv
    ```

## Parameters 

### Action Parameters
Activate at least one action: `--cells`, `--boundaries`, `--summary`, or `--update-catalog` (or `--all`).

- `--cells`: Convert cell points/centroids to PMTiles with counts and top cluster index (topK).
- `--boundaries`: Convert cell polygons to GeoJSON with `cell_id` and `topK` attributes.
- `--summary`: Write a JSON summary of inputs/outputs and basic stats.
- `--update-catalog`: Append generated layers to an existing `catalog.yaml` (not included in `--all`).
- `--all`: Enable `--cells`, `--boundaries`, and `--summary` together.


### Input/Output Parameters
- `--outprefix` (str): Path prefix for outputs (e.g., `/path/to/out/sample`).

=== "Xenium Input"

    - `--in-json` (str): Path to JSON with keys `CELL`, `BOUNDARY`, `CLUSTER`, `DE` (JSON mode).
    - `--indir` (str): Path to Xenium Ranger output directory (manual mode base path).
    - `--csv-cells` (str): Relative path under `--indir` to cells table (CSV/Parquet).
    - `--csv-boundaries` (str): Relative path under `--indir` to cell boundaries (GeoJSON/CSV).
    - `--csv-clust` (str): Relative path under `--indir` to clustering results (CSV).
    - `--csv-diffexp` (str): Relative path under `--indir` to differential expression results (CSV).

=== "Visium HD Input"

    - `--in-json` (str): Path to JSON with keys `CELL_FEATURE_MEX`, `CELL_GEOJSON`, `CLUSTER`, `DE` (JSON mode).
    - `--indir` (str): Path to Space Ranger output directory (manual mode base path).
    - `--mtx-cell` (str): Relative path under `--indir` to cell-feature matrix directory.
    - `--geojson-cells` (str): Relative path under `--indir` to cell geometries (GeoJSON).
    - `--csv-clust` (str): Relative path under `--indir` to clustering results (CSV).
    - `--csv-diffexp` (str): Relative path under `--indir` to differential expression results (CSV).

??? note "Auxiliary Paramaters"

    **PMTiles Conversion**

      - `--min-zoom` (int, default: `10`): Minimum zoom level.
      - `--max-zoom` (int, default: `18`): Maximum zoom level.
      - `--max-tile-bytes` (int, default: `5000000`): Max bytes per tile.
      - `--max-feature-counts` (int, default: `500000`): Max features per tile.
      - `--preserve-point-density-thres` (int, default: `1024`): Density preservation threshold.

    **Coloring and DE**

      - `--tsv-cmap` (str): Color table for clusters (defaults to a fixed color map).
      - `--de-max-pval` (float, default: `0.01`): Max p‑value for DE.
      - `--de-min-fc` (float, default: `1.2`): Min fold change for DE.

    **Catalog Parameters**

      - `--catalog-yaml` (str): Path to catalog to update when using `--update-catalog`. See [Asset Packaging](./run_cartload2.md#parameters) for how layers are organized in the catalog.

    **Environment Parameters**

      - `--tippecanoe` (str): Path to tippecanoe binary (required for PMTiles).
      - `--parquet-tools` (str): Path to parquet-tools if reading Parquet inputs.

    **Run Parameters**

      - `--threads` (int): Threads for tippecanoe conversion.
      - `--log` (flag), `--log-suffix` (str): Write logs alongside outputs.

## Outputs

- Cells PMTiles (`*.cells.pmtiles`): cell centroids/points with attributes `cell_id`, `count`, `topK`.
- Boundaries GeoJSON (`*.boundaries.geojson`): cell polygons with `cell_id`, `topK`.
- Summary TSV/JSON (when `--summary`).
