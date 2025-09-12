# Xenium Pipeline Add-on Modules

## Overview

`CartLoader` provides add-on modules to support the [Xenium pipeline orchestrator](./run_xenium.md):

* `load_xenium_ranger` scans a 10x Xenium Ranger output directory and writes a JSON manifest of detected assets.


## Add-on Modules
=== "Load Xenium Ranger"
    ### Introduction
    `load_xenium_ranger` scans a 10x Xenium Ranger output directory and writes a JSON manifest of detected assets. The manifest consolidates paths for Spatial Gene Expression (SGE), cell analysis outputs (cells, boundaries, clusters, DE, embeddings), and morphology images (OME‑TIFFs).

    Typical next steps use this manifest as input for downstream commands such as `run_xenium`, `sge_convert`, `import_cell`, and `import_image`.

    ---
    ### Requirements

    - Xenium Ranger output directory (unpacked) containing standard files and subfolders.
    - Optional: a writable directory for temporary archive extraction when using `--unzip-dir`.

    ---
    ### Example Usage

    ```bash
    cartloader load_xenium_ranger \
      --in-dir /path/to/xenium/ranger/output \
      --out-json /path/to/out/xenium_ranger_assets.json \
      --unzip-dir /path/to/out/tmp
    ```

    ---
    ### Actions

    Detect standard Xenium Ranger assets under `--in-dir` and build a structured JSON manifest (SGE, CELLS, IMAGES). Paths are resolved from files and supported archives (via `--unzip-dir`); only present items are included:

    - SGE:
        - Required: `TRANSCRIPT` — `transcripts.csv.gz` or `transcripts.parquet`
    - CELLS:
        - Required: `CELL` — `cells.csv.gz` or `cells.parquet`
        - Required: `BOUNDARY` — `cell_boundaries.csv.gz`
        - Optional: `CELL_FEATURE_MEX` — `cell_feature_matrix/` (also `cell_feature_matrix.tar.gz`)
        - Required: `CLUSTER` — `analysis/clustering/gene_expression_graphclust/clusters.csv` (also within `analysis.tar.gz`)
        - Required: `DE` — `analysis/diffexp/gene_expression_graphclust/differential_expression.csv` (also within `analysis.tar.gz`)
        - Optional: `PCA_PROJ` — `analysis/pca/gene_expression_*/projection.csv` (also within `analysis.tar.gz`)
        - Optional: `PCA_VAR` — `analysis/pca/gene_expression_*/variance.csv` (also within `analysis.tar.gz`)
        - Optional: `UMAP_PROJ` — `analysis/umap/gene_expression_*/projection.csv` (also within `analysis.tar.gz`)
    - IMAGES (OME‑TIFF, all optional):
        - `DAPI_OME` — `morphology_focus.ome.tif` or `morphology_focus/morphology_focus_0000.ome.tif`
        - `BOUNDARY_OME` — `morphology_focus/morphology_focus_0001.ome.tif`
        - `INTERIOR_RNA_OME` — `morphology_focus/morphology_focus_0002.ome.tif`
        - `INTERIOR_PROTEIN_OME` — `morphology_focus/morphology_focus_0003.ome.tif`
        - `DAPI_MIP_OME` — `morphology_mip.ome.tif`
      <!-- - Note: `DAPI_3D_OME` may appear in example manifests but is not auto‑detected by this command. -->

    Note: Exact availability varies by dataset and pipeline options. Missing entries are omitted; zip archives listed above are unpacked or scanned when `--unzip-dir` is provided.

    ---
    ### Parameters

    - `--in-dir` (str, required): Xenium Ranger output directory to scan.
    - `--out-json` (str, required): Path to write the generated assets JSON.
    - `--unzip-dir` (str, optional): Directory to extract archives found under `--in-dir` (defaults to `--in-dir`).
    - `--overwrite` (flag): Overwrite `--out-json` if it already exists.

    ---
    ### Output

    Writes a JSON manifest at `--out-json`. 
    
    Example: [`includes/xenium_ranger_assets.human_lung_cancer.json`](../../includes/xenium_ranger_assets.human_lung_cancer.json)

    ```json
    {% include-markdown "../../includes/xenium_ranger_assets.human_lung_cancer.json" %}
    ```
