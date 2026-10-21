# Xenium Pipeline Add-on Modules (Legacy)

!!! warning "Legacy"
    `load_xenium_ranger` supports the legacy [`run_xenium`](./run_xenium.md) orchestrator, which is superseded by [`run_together`](./run_together.md); `run_together` detects Xenium inputs from `--in-dir` through its [Xenium profile](./platforms/xenium.md). This page is kept for reference and is not part of the recommended workflow; see [Legacy Orchestrators](./legacy.md).

## Overview

`load_xenium_ranger` scans a 10x Xenium Ranger output directory and writes a JSON manifest of detected assets (SGE, cell analysis outputs, and morphology images). The manifest can be passed as `--in-json` to `sge_convert`, [`import_xenium_cell`](./import_cell.md), and [`import_image`](./import_image.md), and is read by `run_xenium`. It is also a way to unpack the Xenium `analysis.tar.gz` and `cell_feature_matrix.tar.gz` archives, which `run_together` does not extract.

---
## Requirements

- A Xenium Ranger output directory.
- Optional: a writable directory for archive extraction (`--unzip-dir`).

---
## Actions

### Load Xenium Ranger (`load_xenium_ranger`)

For each key below, the first existing path under `--in-dir` is used. If none exists and the key lists an archive, the archive (under `--in-dir`) is extracted into `--unzip-dir` and searched there.

- SGE:
    - Required: `TRANSCRIPT` — `transcripts.csv.gz` or `transcripts.parquet`
- CELLS:
    - Required: `CELL` — `cells.csv.gz` or `cells.parquet`
    - Required: `BOUNDARY` — `cell_boundaries.csv.gz`
    - Optional: `CELL_FEATURE_MEX` — `cell_feature_matrix/` (archive: `cell_feature_matrix.tar.gz`)
    - Required: `CLUSTER` — `analysis/clustering/gene_expression_graphclust/clusters.csv` (archive: `analysis.tar.gz`)
    - Required: `DE` — `analysis/diffexp/gene_expression_graphclust/differential_expression.csv` (archive: `analysis.tar.gz`)
    - Optional: `PCA_PROJ` — `analysis/pca/gene_expression_10_components/projection.csv` (archive: `analysis.tar.gz`)
    - Optional: `PCA_VAR` — `analysis/pca/gene_expression_10_components/variance.csv` (archive: `analysis.tar.gz`)
    - Optional: `UMAP_PROJ` — `analysis/umap/gene_expression_2_components/projection.csv` (archive: `analysis.tar.gz`)
- IMAGES (OME‑TIFF, all optional):
    - `DAPI_OME` — `morphology_focus.ome.tif` or `morphology_focus/morphology_focus_0000.ome.tif`
    - `BOUNDARY_OME` — `morphology_focus/morphology_focus_0001.ome.tif`
    - `INTERIOR_RNA_OME` — `morphology_focus/morphology_focus_0002.ome.tif`
    - `INTERIOR_PROTEIN_OME` — `morphology_focus/morphology_focus_0003.ome.tif`
    - `DAPI_MIP_OME` — `morphology_mip.ome.tif` (present in older Xenium outputs only)

The 3D DAPI stack (`morphology.ome.tif`) is not detected.

!!! note "Missing files"
    Paths are exact; availability varies by dataset. Missing optional entries are left out of the manifest, but a missing required entry stops the command with an error listing what was not found. An optional entry whose archive is listed above also stops the command if neither the file nor the archive exists (e.g., `cell_feature_matrix/` and `cell_feature_matrix.tar.gz` both absent, or a missing `UMAP_PROJ` without `analysis.tar.gz`).

If `--out-json` already exists and `--overwrite` is not set, its entries are kept and only missing keys are resolved.

---
## Example Usage

```bash
cartloader load_xenium_ranger \
  --in-dir /path/to/xenium/ranger/output \
  --out-json /path/to/out/xenium_ranger_assets.json \
  --unzip-dir /path/to/out/unzip
```

---
## Parameters

- `--in-dir` (str, required): Xenium Ranger output directory to scan.
- `--out-json` (str, required): Path to write the assets JSON.
- `--unzip-dir` (str): Directory to extract archives into (default: `--in-dir`).
- `--overwrite` (flag): Re-resolve all entries instead of keeping those in an existing `--out-json`.

---
## Output

Writes a JSON manifest at `--out-json`.

Example (`includes/xenium_ranger_assets.human_lung_cancer.json`):

```json
{% include-markdown "../../includes/xenium_ranger_assets.human_lung_cancer.json" %}
```
