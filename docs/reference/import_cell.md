# Cell Analysis Import

## Overview

Import platform‑specific cell analysis outputs (cells, boundaries, clusters/DE), convert them to web‑ready layers (PMTiles/GeoJSON), and optionally register them in a `CartLoader` catalog. Supported platforms include 10x Xenium and 10x Visium HD.

---
## Requirements

- Cell analysis outputs from the platform (see per‑platform sections).
- For PMTiles conversion: `tippecanoe` installed and on PATH.
- Optional catalog update: an existing `catalog.yaml` from [Asset Packaging](./run_cartload2.md).

---
## Actions

!!! warning "Action Specifications"
    No action runs by default. Activate at least one using the [action parameters](#action-parameters).

### Cell Conversion Step (`--cells`)
Convert cell points/centroids, clusters, and differential expression results to PMTiles with counts and top cluster index (topK). A JSON summary of inputs/outputs and basic stats will be written.

### Cell Boundary Conversion Step (`--boundaries`)
Convert cell polygons to GeoJSON with `cell_id` and `topK` attributes. A JSON summary of inputs/outputs and basic stats will be written.

### Catalog Update (`--update-catalog`)
Append generated layers to an existing `catalog.yaml`.

### UMAP Conversion Step (`--umap`)
Convert UMAP embeddings to PMTiles.

--- 
## Example Usage
!!! warning "Two Input Modes"

    `CartLoader` accepts two input modes for each platform:

    - JSON mode: pass `--in-json` with required keys.
    - Manual mode: pass `--in-dir` with specific file locations.

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
      --in-dir /path/to/xenium_ranger_output \
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
      --cells \
      --boundaries \
      --outprefix /path/to/out/spaceranger \
      --in-json /path/to/visiumhd_cells.json

    # Manual input mode
    cartloader import_visiumhd_cell \
      --cells \
      --boundaries \
      --outprefix /path/to/out/spaceranger \
      --in-dir /path/to/spaceranger_output \
      --mtx-cell filtered_feature_bc_matrix \
      --geojson-cells spatial/cells.geojson \
      --csv-clust analysis/clustering/gene_expression_graphclust/clusters.csv \
      --csv-diffexp analysis/diffexp/gene_expression_graphclust/differential_expression.csv
    ```

---

## Parameters

### Action Parameters

- `--cells` (flag): Import [segmented cells](#cell-conversion-step---cells) and generate PMTiles.
- `--boundaries` (flag): Import [segmented cell boundaries](#cell-boundary-conversion-step---boundaries) and generate GeoJSON/PMTiles.
- `--umap` (flag): Import [UMAP projection](#umap-conversion-step---umap) and generate PMTiles.
- `--update-catalog` (flag): [Update an existing catalog.yaml](#catalog-update---update-catalog).
- `--all` (flag): Enable all actions (`--cells`, `--boundaries`, `--umap`).

### Input/Output Parameters
- `--outprefix` (str, required): Output prefix (e.g., `/path/to/out/sample`).
- `--id` (str): Identifier for the cell factor (default: basename of `--outprefix`).
- `--name` (str): Display name for the cell factor (default: basename of `--outprefix`).
- `--tmp-dir` (str): Temporary directory for intermediate files.

=== "Xenium Input"
    **JSON Mode**:
    - `--in-json` (str): Path to JSON with keys `CELL`, `BOUNDARY`, `CLUSTER`, `DE`, `UMAP_PROJ`.

    **Manual Mode** (relative to `--in-dir`):
    - `--in-dir` (str): Base input directory.
    - `--csv-cells` (str): Cell table (CSV/Parquet) (default: `cells.csv.gz`).
    - `--csv-boundaries` (str): Cell boundaries (CSV) (default: `cell_boundaries.csv.gz`).
    - `--csv-clust` (str): Clustering results (CSV) (default: `analysis/clustering/gene_expression_graphclust/clusters.csv`).
    - `--csv-diffexp` (str): Differential expression (CSV) (default: `analysis/diffexp/gene_expression_graphclust/differential_expression.csv`).
    - `--csv-umap` (str): UMAP projection (CSV) (default: `analysis/pca/gene_expression_10_components/projection.csv`).

=== "Visium HD Input"
    **JSON Mode**:
    - `--in-json` (str): Path to JSON with keys `CELL_FEATURE_MEX`, `CELL_GEOJSON`, `CLUSTER`, `DE`, `UMAP_PROJ`.

    **Manual Mode** (relative to `--in-dir`):
    - `--in-dir` (str): Base input directory.
    - `--mtx-cells` (str): Cell-feature matrix directory (default: `segmented_outputs/filtered_feature_cell_matrix`).
    - `--geojson-cells` (str): Cell segmentations (GeoJSON) (default: `segmented_outputs/cell_segmentations.geojson`).
    - `--csv-clust` (str): Clustering results (CSV) (default: `analysis/clustering/gene_expression_graphclust/clusters.csv`).
    - `--csv-diffexp` (str): Differential expression (CSV) (default: `analysis/diffexp/gene_expression_graphclust/differential_expression.csv`).
    - `--csv-umap` (str): UMAP projection (CSV) (default: `analysis/pca/gene_expression_10_components/projection.csv`).
    - `--scale-json` (str): Scale factors JSON (default: `spatial/scalefactors_json.json`, or None).
    - `--units-per-um` (float): Coordinate units per µm (default: 1).

??? note "Auxiliary Parameters"

    **PMTiles Conversion**:
      - `--min-zoom` (int, default: 10): Min zoom for cells/boundaries.
      - `--max-zoom` (int, default: 18): Max zoom for cells/boundaries.
      - `--umap-min-zoom` (int, default: 0): Min zoom for UMAP.
      - `--umap-max-zoom` (int, default: 18): Max zoom for UMAP.
      - `--max-tile-bytes` (int, default: 5000000): Max bytes per tile.
      - `--max-feature-counts` (int, default: 500000): Max features per tile.
      - `--preserve-point-density-thres` (int, default: 1024): Point density threshold.

    **Coloring and DE**:
      - `--tsv-cmap` (str): Color table TSV (default: `assets/fixed_color_map_60.tsv`).
      - `--de-max-pval` (float, default: 0.01): Max p-value.
      - `--de-min-fc` (float, default: 1.2): Min fold change.

    **Catalog**:
      - `--catalog-yaml` (str): Input catalog to update.
      - `--out-catalog-yaml` (str): Output catalog path (if different from input).

    **Column Names (Collapsible)**:
      *Xenium Only*:
      - `--cells-colname-cell-id` (str, default: `cell_id`)
      - `--cells-colname-x` (str, default: `x_centroid`)
      - `--cells-colname-y` (str, default: `y_centroid`)
      - `--cells-colname-count` (str, default: `transcript_counts`)
      - `--boundaries-colname-cell-id` (str, default: `cell_id`)
      - `--boundaries-colname-x` (str, default: `vertex_x`)
      - `--boundaries-colname-y` (str, default: `vertex_y`)
      
      *Common/Both*:
      - `--clust-colname-barcode` (str, default: `Barcode`)
      - `--clust-colname-cluster` (str, default: `Cluster`)
      - `--umap-colname-barcode` (str, default: `Barcode`)
      - `--umap-colname-x` (str, default: `UMAP-1`)
      - `--umap-colname-y` (str, default: `UMAP-2`)

    **Environment**:
      - `--tippecanoe` (str): Path to tippecanoe binary.
      - `--parquet-tools` (str): Path to parquet-tools (Xenium only).
      - `--R` (str): Path to R/Rscript binary.
      - `--threads` (int): Number of threads.
      - `--log` (flag): Enable logging.
      - `--log-suffix` (str): Log file suffix.

---
## Output

- Cells PMTiles (`*.cells.pmtiles`): cell centroids/points with attributes `cell_id`, `count`, `topK`.
- Boundaries GeoJSON (`*.boundaries.geojson`): cell polygons with `cell_id`, `topK`.
- Summary TSV/JSON (when `--summary`).
