# Xenium Pipeline

## Overview

An end‑to‑end workflow for 10x Xenium data using CartLoader: convert inputs, optionally reorient/stitch, run FICTURE, package assets, and (optionally) upload for sharing.

## Input Data

This tutorial uses SGE data generated with the 10x Genomics Xenium platform.

**SGE**

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_xenium.md"
%}

## Inputs

=== "JSON Mode"

    Provide a JSON file describing Xenium assets (cells, boundaries, clusters, DE).

    ```json
    {
      "CELL": "/path/to/xenium/analysis/cell_feature_matrix/cells.csv",
      "BOUNDARY": "/path/to/xenium/analysis/cell_segmentation/cells_boundary.geojson",
      "CLUSTER": "/path/to/xenium/analysis/clustering/gene_expression_graphclust/clusters.csv",
      "DE": "/path/to/xenium/analysis/diffexp/gene_expression_graphclust/differential_expression.csv"
    }
    ```

=== "Manual Mode"

    Use an output directory and platform‑specific relative paths.

    ```bash
    # Example relative paths under --indir
    analysis/cell_feature_matrix/cells.csv
    analysis/cell_segmentation/cells_boundary.geojson
    analysis/clustering/gene_expression_graphclust/clusters.csv
    analysis/diffexp/gene_expression_graphclust/differential_expression.csv
    ```

## Steps (copy/paste examples)

- sge_convert: Prepare FICTURE‑ready SGE
- run_ficture2: Train and decode
- run_cartload2: Package assets (PMTiles + catalog)
- import_xenium_cell: Add cell layers (points + boundaries)
- upload (optional): AWS or Zenodo

```bash
# Cell layers (JSON mode)
cartloader import_xenium_cell \
  --cells --boundaries --summary \
  --outprefix /path/to/out/xenium \
  --in-json /path/to/xenium_cells.json

# Packaging
cartloader run_cartload2 \
  --fic-dir /path/to/run_ficture2/results \
  --out-dir /path/to/cartload2 \
  --id xenium_example
```

## Outputs

- PMTiles for molecules/factors/rasters
- Cell point PMTiles and boundary GeoJSON
- Catalog YAML with all layers

## See Also

- Reference: import_cell.md
- Reference: run_ficture2.md
- Reference: run_cartload2.md
