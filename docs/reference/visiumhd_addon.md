# Visium HD Pipeline Add-on Modules (Legacy)

!!! warning "Legacy"
    `load_space_ranger` supports the legacy [`run_visiumhd`](./run_visiumhd.md) orchestrator, which is superseded by [`run_together`](./run_together.md); `run_together` detects Visium HD inputs from `--in-dir` through its [Visium HD profile](./platforms/visium_hd.md). This page is kept for reference and is not part of the recommended workflow; see [Legacy Orchestrators](./legacy.md).

## Overview

`load_space_ranger` scans a 10x Visium HD Space Ranger output directory and writes a JSON manifest of detected assets: SGE inputs (2 µm bins and spatial metadata), segmented cell analysis outputs, the 8 µm and 16 µm square-bin analyses, and an H&E image. The manifest can be passed as `--in-json` to `sge_convert`, [`import_visiumhd_cell`](./import_cell.md), [`import_visiumhd_square`](./import_square.md), and [`import_image`](./import_image.md), and is read by `run_visiumhd`. It is also a way to unpack Space Ranger `.tar.gz` archives, which `run_together` does not extract.

---
## Requirements

- A Space Ranger output directory.
- Optional: a writable directory for archive extraction (`--unzip-dir`).

---
## Actions

### Load Space Ranger (`load_space_ranger`)

For each key below, the first existing path under `--in-dir` is used. If none exists, any matching archive found anywhere under `--in-dir` (`*_segmented_outputs.tar.gz` for cells; `*_square_<NNN>um_binned_outputs.tar.gz` or `*_binned_outputs.tar.gz` for bins) is extracted into `--unzip-dir` and searched there.

- SGE (required keys):
    - `TRANSCRIPT_MEX`: `square_002um/filtered_feature_bc_matrix` or `binned_outputs/square_002um/filtered_feature_bc_matrix`
    - `SCALE`: `square_002um/spatial/scalefactors_json.json` or `binned_outputs/square_002um/spatial/scalefactors_json.json`
    - `POSITION`: `square_002um/spatial/tissue_positions.parquet` or `binned_outputs/square_002um/spatial/tissue_positions.parquet`

- CELLS (segmented outputs):
    - Required: `CELL_FEATURE_MEX` — `segmented_outputs/filtered_feature_cell_matrix`
    - Required: `CELL_GEOJSON` — `segmented_outputs/cell_segmentations.geojson`
    - Optional: `NUCLEUS_GEOJSON` — `segmented_outputs/nucleus_segmentations.geojson`
    - Required: `CLUSTER` — `segmented_outputs/analysis/clustering/gene_expression_graphclust/clusters.csv`
    - Required: `DE` — `segmented_outputs/analysis/diffexp/gene_expression_graphclust/differential_expression.csv`
    - Optional: `PCA_PROJ` — `segmented_outputs/analysis/pca/gene_expression_10_components/projection.csv`
    - Optional: `PCA_VAR` — `segmented_outputs/analysis/pca/gene_expression_10_components/variance.csv`
    - Optional: `UMAP_PROJ` — `segmented_outputs/analysis/umap/gene_expression_2_components/projection.csv`

- GRID_8um and GRID_16um (square bins of 8 µm and 16 µm; all optional), each under `square_008um/` or `binned_outputs/square_008um/` (`square_016um/` for 16 µm):
    - `GRID_FEATURE_MEX`: `filtered_feature_bc_matrix`
    - `SCALE`: `spatial/scalefactors_json.json`
    - `POSITION`: `spatial/tissue_positions.parquet`
    - `CLUSTER`: `analysis/clustering/gene_expression_graphclust/clusters.csv`
    - `DE`: `analysis/diffexp/gene_expression_graphclust/differential_expression.csv`
    - `PCA_PROJ`: `analysis/pca/gene_expression_10_components/projection.csv`
    - `PCA_VAR`: `analysis/pca/gene_expression_10_components/variance.csv`
    - `UMAP_PROJ`: `analysis/umap/gene_expression_2_components/projection.csv`

- IMAGES (optional):
    - `HnE`: The first file under `--in-dir` (searched recursively) ending in `_tissue_image.btf`, `_tissue_image.tif`, or `_tissue_image.tiff`, in that order of preference. Archives are not searched for images.

!!! note "Missing files"
    Paths are exact; availability varies by dataset. Missing optional entries are left out of the manifest (a section with nothing found is written empty), but a missing required entry stops the command with an error listing what was not found.

If `--out-json` already exists and `--overwrite` is not set, its entries are kept and only missing keys are resolved.

---
## Example Usage

```bash
cartloader load_space_ranger \
  --in-dir /path/to/SpaceRanger_*_outs \
  --out-json /path/to/out/visiumhd_space_ranger_assets.json
```

---
## Parameters

- `--in-dir` (str, required): Space Ranger output directory to scan.
- `--out-json` (str, required): Path to write the assets JSON.
- `--unzip-dir` (str): Directory to extract archives into (default: `--in-dir`).
- `--overwrite` (flag): Re-resolve all entries instead of keeping those in an existing `--out-json`.

---
## Output

Writes a JSON manifest at `--out-json`.

Example (`includes/visiumhd_space_ranger_assets.visiumhd_3prime_mouse_brain.json`):

```json
{% include-markdown "../../includes/visiumhd_space_ranger_assets.visiumhd_3prime_mouse_brain.json" %}
```
