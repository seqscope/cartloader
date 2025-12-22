# Visium HD Pipeline Orchestrator


## Overview

An end‑to‑end workflow for 10x Visium HD data using `CartLoader` to run steps (including load Space Ranger assets, convert SGE, run FICTURE, import cells, square-bin outputs with optional UMAPs, background images, package assets, and upload for sharing) as needed.

---
## Requirements

- Space Ranger output
- Pre-installed tools depending on the actions: `spatula`, `punkst`, `gzip`, `sort`, `python`, `go-pmtiles`, `gdal`, `tippecanoe`, `parquet-tools`, `aws`

---
## Example Usage

=== "Run all modules together"

    ```bash
    cartloader run_visiumhd \
      --load-space-ranger \
      --sge-convert \
      --run-ficture2 \
      --import-cells \
      --import-images \
      --run-cartload2 \
      --upload-aws \
      --space-ranger-dir /path/to/space/ranger/output \
      --out-dir /path/to/out/dir \
      --s3-dir s3://example_bucket/example-id \
      --width 12 \
      --n-factor 12 \
      --id example-id \
      --spatula /path/to/spatula/binary \
      --ficture2 /path/to/punkst/directory \  
      --pmtiles /path/to/pmtiles/binary \
      --tippecanoe /path/to/tippecanoe/binary \
      --aws /path/to/aws/cli/binary \
      --n-jobs 10 \
      --threads 10
    ```

---
## Actions

!!! info "Action Specifications"
    
    `run_visiumhd` orchestrates multiple `CartLoader` modules; enable and combine actions with flags.
    
    No action runs by default. Activate at least one action flag.

See details for each flag:

- [`--load-space-ranger`](./load_space_ranger.md#actions)
- [`--sge-convert`](./sge_convert.md#actions)
- [`--run-ficture2`](./run_ficture2.md#actions)
- [`--import-cells`](./import_cell.md#actions)
- `--import-squares`: Import Space Ranger square bins (multi-resolution) including optional UMAP coordinates per bin size; writes square assets JSONs and updates `catalog.yaml`.
- [`--import-images`](./import_image.md#actions)
- [`--run-cartload2`](./run_cartload2.md#actions)
- [`--upload-aws`](./upload_aws.md#actions)
- [`--upload-zenodo`](./upload_zenodo.md#actions)

---
## Parameters

Below are the core arguments you’ll typically set. For all other options, expand the collapsible "Auxiliary Parameters" section.

### Action Parameters

See [Action](#action).

### Input/Output Parameters

- `--space-ranger-dir`: Space Ranger output directory (required for manual mode or with `--load-space-ranger`)
- `--out-dir`: Output root directory; writes `sge/`, `ficture2/`, `cartload2/` (required)
- `--space-ranger-assets`: Path to assets JSON from `--load-space-ranger` (default: `<out_dir>/space_ranger_assets.json`)
- `--ext-fic-dir`: External FICTURE directory to import when using `--import-ext-ficture2`

### SGE Convert (with `--sge-convert`)

- `--units-per-um`: Coordinate units per micron in raw SGE (float)
- `--exclude-feature-regex`: Regex for features to exclude (e.g., controls)
- `--filter-by-density`: Enable density-based filtering (bool)

If using JSON mode (with `--load-space-ranger`), SGE inputs are discovered from `--space-ranger-assets`.
In manual mode, provide inputs relative to `--space-ranger-dir`:
- `--mex-transcript`: Transcript MEX directory (`square_002um/filtered_feature_bc_matrix`)
- `--parquet-position`: Path to `square_002um/spatial/tissue_positions.parquet`
- `--json-scale` or `--units-per-um`: Pixel‑to‑µm scale source

### FICTURE (with `--run-ficture2`)
- `--width`: Comma‑separated hexagon widths in µm for train/projection (required)
- `--n-factor`: Comma‑separated factor counts for training (required)
- `--colname-feature`: Feature column name (default: `gene`)
- `--colname-count`: Count column name (default: `count`)
- `--fic-include-feature-regex` (str): Regex of feature names to include in FICTURE analysis.
- `--fic-exclude-feature-regex` (str): Regex of feature names to exclude in FICTURE analysis.

### Import Cells (with `--import-cells`)

- `--cell-id`: Asset ID/prefix for cell outputs
- `--cell-name`: Human‑readable name (defaults to `--cell-id`)
- `--tsv-cmap`: TSV colormap for cluster colors

Manual mode inputs (relative to `--space-ranger-dir`):
- `--csv-clust`: Cluster assignments CSV
- `--csv-diffexp`: Differential expression CSV

### Import Squares (with `--import-squares`)

- `--square-id`: Prefix for square-bin assets (default: `spaceranger`)
- `--use-parquet-tools`: Use `parquet-tools` instead of polars/pigz for Parquet→CSV conversion (slower on large files)
- `--square-input`: Manual mode only: one or more square bin definitions, each as `bin_size,bin_pos_parquet,bin_csv_cluster,bin_csv_diffexp,bin_scale_json,bin_csv_umap`. UMAP CSV is optional; include to expose 2D embeddings per bin size.

### Import Images (with `--import-images`)

- `--image-ids`: One or more image IDs to import (e.g., `HnE_BTF`)
- `--all-images`: Import all images listed in `--space-ranger-assets`
- `--transparent-below`: Make pixels below value transparent (0–255; default: 1)

Manual mode inputs (relative to `--space-ranger-dir`):
- `--tifs`: One or more BTF/TIFF paths

### CartLoad2 (with `--run-cartload2`)

- `--id`: Catalog ID; no whitespace; prefer `-` over `_` (required)
- `--title`: Human‑readable title (quote if spaces)
- `--desc`: Short description (quote if spaces)

### Upload (with `--upload-aws` / `--upload-zenodo`)

- `--s3-dir`: S3 destination, e.g., `s3://bucket/prefix` (required for AWS)
- `--zenodo-token`: Path to file containing Zenodo token (required for Zenodo)
- `--zenodo-deposition-id`: Existing deposition ID (omit to create new)
- `--zenodo-title`: Deposition title (defaults to `--title`)
- `--creators`: One or more `"Lastname, Firstname"` entries (each quoted)

??? info "Auxiliary Parameters"

    #### Run Options

    - `--dry-run`: Print commands without executing (bool; default: false)
    - `--restart`: Ignore existing outputs and rerun steps (bool; default: false)
    - `--n-jobs`, `-j`: Parallel jobs to run (int; default: 1)
    - `--threads`: Max threads per job, e.g., tippecanoe (int; default: 4)

    #### Environment / Tools

    - `--pmtiles`: Path to `pmtiles` binary (default: `pmtiles`)
    - `--gdal_translate`: Path to `gdal_translate` (default: `gdal_translate`)
    - `--gdaladdo`: Path to `gdaladdo` (default: `gdaladdo`)
    - `--gdalwarp`: Path to `gdalwarp` (default: `gdalwarp`)
    - `--tippecanoe`: Path to `tippecanoe` (default: repo submodule)
    - `--spatula`: Path to `spatula` (default: `spatula`)
    - `--parquet-tools`: Path to `parquet-tools` (default: `parquet-tools`)
    - `--ficture2`: Path to punkst/FICTURE repo (default: repo submodule)
    - `--aws`: Path to AWS CLI (default: `aws`)

---
## Output

Outputs depend on the enabled action flags. Refer to each module page for full output details and file formats.

- `--load-space-ranger`: A Space Ranger assets JSON summarizing detected inputs (path: `--space-ranger-assets`). See an example in the Visium HD pipeline vignette.
- `--sge-convert`: Unified SGE files and an assets manifest (`sge_assets.json`). See details in [sge_convert reference](./sge_convert.md#output)
- `--run-ficture2` or `--import-ext-ficture2`: Generates FICTURE analysis artifacts (factor models, decoded maps, PMTiles, summaries). See details in [run_ficture2 reference](./run_ficture2.md#output)
- `--import-cells`: Cells PMTiles, boundaries GeoJSON, and summaries under `cartload2/`. See details in [import_cells reference](./import_cell.md#outputs).
- `--import-squares`: Square-bin PMTiles plus assets JSONs per bin size (e.g., `<square-id>-sq050_assets.json`), including optional UMAP TSVs when provided.
- `--import-images`: Background image PMTiles (and intermediates) under `cartload2/`. See details in [import_images reference](./import_image.md#output).
- `--run-cartload2`: Sources packaged into PMTiles and a `catalog.yaml`. See details in [run_cartload2 reference](./run_cartload2.md#output).
- `--upload-aws`: Uploads generated assets and `catalog.yaml` to the specified S3 path.
- `--upload-zenodo`: Uploads `catalog.yaml` and selected assets to Zenodo.
