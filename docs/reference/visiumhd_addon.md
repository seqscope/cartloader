# Visium HD Pipeline Add-ons Modules

## Overview

`CartLoader` provides add-on modules to support [Visium HD pipeline orchestrator](./run_visiumhd.md):

* `load_space_ranger` scans a 10x Visium HD Space Ranger output directory and writes a JSON manifest of detected assets.

## Add-on Modules
=== "Load Space Ranger"
    ### Introduction

    Create a manifest consolidates locations for SGE inputs (binned features and spatial metadata), segmented cell analysis outputs, hex-binned analysis at 8 µm and 16 µm, and optional tissue images. 

    Typical next steps use this manifest as input for `sge_convert`, `import_cell`, `import_image`, or an end‑to‑end pipeline.
    
    ---

    ### Requirements

    - Space Ranger output directory (unpacked) containing standard files and subfolders.
    - Optional: a writable directory for temporary archive extraction when using `--unzip-dir`.
    
    ---
    ### Example Usage

    ```bash
    cartloader load_space_ranger \
      --in-dir /path/to/SpaceRanger_*_outs \
      --out-json /path/to/out/visiumhd_space_ranger_assets.json 
    ```

    ---
    ### Actions

    Detect standard Space Ranger assets under `--in-dir` and build a structured JSON manifest. Paths are resolved from files and supported archives (via `--unzip-dir`); only present items are included.

    - SGE (required keys):
        - `TRANSCRIPT_MEX`: `square_002um/filtered_feature_bc_matrix` or `binned_outputs/square_002um/filtered_feature_bc_matrix` (also in `*_square_002um_binned_outputs.tar.gz` or `*_binned_outputs.tar.gz`)
        - `SCALE`: `square_002um/spatial/scalefactors_json.json` or `binned_outputs/square_002um/spatial/scalefactors_json.json` (same archives)
        - `POSITION`: `square_002um/spatial/tissue_positions.parquet` or `binned_outputs/square_002um/spatial/tissue_positions.parquet` (same archives)

    - CELLS (segmented outputs):
        - Required: `CELL_FEATURE_MEX` — `segmented_outputs/filtered_feature_cell_matrix` (also in `*_segmented_outputs.tar.gz`)
        - Required: `CELL_GEOJSON` — `segmented_outputs/cell_segmentations.geojson` (also in `*_segmented_outputs.tar.gz`)
        - Optional: `NUCLEUS_GEOJSON` — `segmented_outputs/nucleus_segmentations.geojson` (also in `*_segmented_outputs.tar.gz`)
        - Required: `CLUSTER` — `segmented_outputs/analysis/clustering/gene_expression_graphclust/clusters.csv` (also in `*_segmented_outputs.tar.gz`)
        - Required: `DE` — `segmented_outputs/analysis/diffexp/gene_expression_graphclust/differential_expression.csv` (also in `*_segmented_outputs.tar.gz`)
        - Optional: `PCA_PROJ` — `segmented_outputs/analysis/pca/gene_expression_10_components/projection.csv` (also in `*_segmented_outputs.tar.gz`)
        - Optional: `PCA_VAR` — `segmented_outputs/analysis/pca/gene_expression_10_components/variance.csv` (also in `*_segmented_outputs.tar.gz`)
        - Optional: `UMAP_PROJ` — `segmented_outputs/analysis/umap/gene_expression_2_components/projection.csv` (also in `*_segmented_outputs.tar.gz`)

    - GRID_8um and GRID_16um (hex‑binned; all optional):
        - `GRID_FEATURE_MEX`: `<square_008um|square_016um>/filtered_feature_bc_matrix` or `binned_outputs/<square_...>/filtered_feature_bc_matrix` (also in `*_square_00[8|16]um_binned_outputs.tar.gz` or `*_binned_outputs.tar.gz`)
        - `SCALE`: `<square_...>/spatial/scalefactors_json.json` (same archive support)
        - `POSITION`: `<square_...>/spatial/tissue_positions.parquet` (same archive support)
        - `CLUSTER`: `<square_...>/analysis/clustering/gene_expression_graphclust/clusters.csv` (same archive support)
        - `DE`: `<square_...>/analysis/diffexp/gene_expression_graphclust/differential_expression.csv` (same archive support)
        - `PCA_PROJ`: `<square_...>/analysis/pca/gene_expression_10_components/projection.csv` (same archive support)
        - `PCA_VAR`: `<square_...>/analysis/pca/gene_expression_10_components/variance.csv` (same archive support)
        - `UMAP_PROJ`: `<square_...>/analysis/umap/gene_expression_2_components/projection.csv` (same archive support)

    - IMAGES (optional):
        - `HnE_BTF`: Matches files with suffix `_tissue_image.btf` (also found within `*_square_002um_binned_outputs.tar.gz` or `*_binned_outputs.tar.gz`)

    Note: Exact availability varies by dataset and pipeline options. Missing entries are omitted.

    ---
    ### Parameters

    - `--in-dir` (str, required): Space Ranger output directory to scan.
    - `--out-json` (str, required): Path to write the generated assets JSON.
    - `--unzip-dir` (str, optional): Directory to extract archives found under `--in-dir` (defaults to `--in-dir`).
    - `--overwrite` (flag): Overwrite `--out-json` if it already exists.

    ---
    ## Output

    Writes a JSON manifest at `--out-json`.
    
    Example: [`includes/visiumhd_space_ranger_assets.visiumhd_3prime_mouse_brain.json`](../../includes/visiumhd_space_ranger_assets.visiumhd_3prime_mouse_brain.json)

    ```json
    {% include-markdown "../../includes/visiumhd_space_ranger_assets.visiumhd_3prime_mouse_brain.json" %}
    ```
