# Visium HD Pipeline

## Overview

An end‑to‑end workflow for 10x Visium HD data using CartLoader: convert inputs, optionally reorient/stitch, run FICTURE, package assets, and (optionally) upload for sharing.

## Inputs

=== "JSON Mode"

    Provide a JSON file describing Visium HD assets (cell feature matrix, cell geometries, clusters, DE).

    ```json
    {
      "CELL_FEATURE_MEX": "/path/to/spaceranger/filtered_feature_bc_matrix",
      "CELL_GEOJSON": "/path/to/spaceranger/spatial/cells.geojson",
      "CLUSTER": "/path/to/spaceranger/analysis/clustering/gene_expression_graphclust/clusters.csv",
      "DE": "/path/to/spaceranger/analysis/diffexp/gene_expression_graphclust/differential_expression.csv"
    }
    ```

=== "Manual Mode"

    Use an output directory and platform‑specific relative paths.

    ```bash
    # Example relative paths under --indir
    filtered_feature_bc_matrix
    spatial/cells.geojson
    analysis/clustering/gene_expression_graphclust/clusters.csv
    analysis/diffexp/gene_expression_graphclust/differential_expression.csv
    ```

## Steps (copy/paste examples)

- sge_convert: Prepare FICTURE‑ready SGE
- run_ficture2: Train and decode
- run_cartload2: Package assets (PMTiles + catalog)
- import_visiumhd_cell: Add cell layers (points + boundaries)
- upload (optional): AWS or Zenodo

```bash
# Cell layers (JSON mode)
cartloader import_visiumhd_cell \
  --cells --boundaries --summary \
  --outprefix /path/to/out/visiumhd \
  --in-json /path/to/visiumhd_cells.json

# Packaging
cartloader run_cartload2 \
  --fic-dir /path/to/run_ficture2/results \
  --out-dir /path/to/cartload2 \
  --id visiumhd_example
```

## Outputs

- PMTiles for molecules/factors/rasters
- Cell point PMTiles and boundary GeoJSON
- Catalog YAML with all layers

## See Also

- Reference: import_cell.md
- Reference: run_ficture2.md
- Reference: run_cartload2.md
