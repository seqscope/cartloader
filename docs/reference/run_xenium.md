# Xenium Pipeline Orchestrator

## Overview

An end‑to‑end workflow for 10x Xenium data using `CartLoader` to run steps (including convert inputs, run FICTURE, convert cell analysis results, convert background images, package assets, and upload for sharing) as needed.

---
## Requirements

- Xenium Ranger output
- Pre-installed tools depending on the actions: `spatula`, `punkst`, `gzip`, `sort`, `python`, `go-pmtiles`, `gdal`, `tippecanoe`, `parquet-tools`, `aws`

---
## Actions

!!! warning "Action Specifications"
    `run_xenium` orchestrates multiple `CartLoader` modules; enable and combine actions with flags.
    
    No action runs by default. Activate at least one action flag.

Click to see action details per flag:

* [`--load-xenium-ranger`](./load_xenium_ranger.md#actions)
* [`--sge-convert`](./sge_convert.md#actions)
* [`--run-ficture2`](./run_ficture2.md#actions)
* [`--import-cells`](./import_cell.md#actions)
* [`--import-images`](./import_image.md#actions)
* [`--run-cartload2`](./run_cartload2.md#actions)
* [`--upload-aws`](./upload_aws.md#actions)
* [`--upload-zenodo`](./upload_zenodo.md#actions)

---
## Example Usage

!!! warning "Replace placeholders"

    Replace placeholder paths (e.g., `/path/to/xenium/ranger/output`) or example values (e.g., `example-id` and `12`) before execution.

=== "All modules (end‑to‑end)"

    The example below imports only the `OME_DAPI` image for illustration.

    ```bash
    cartloader run_xenium \
      --load-xenium-ranger \
      --sge-convert \
      --run-ficture2 \
      --import-cells \
      --import-images \
      --run-cartload2 \
      --upload-aws \
      --xenium-ranger-dir /path/to/xenium/ranger/output \
      --out-dir /path/to/out/dir \
      --s3-dir s3://example_bucket/example-id \
      --width 12 \
      --n-factor 12 \
      --id example-id \
      --image-ids OME_DAPI \
      --spatula /path/to/spatula/binary \
      --ficture2 /path/to/punkst/directory \
      --pmtiles /path/to/pmtiles/binary \
      --tippecanoe /path/to/tippecanoe/binary \
      --aws /path/to/aws/cli/binary \
      --n-jobs 10 \
      --threads 10
    ```

=== "FICTURE training + packaging"

    ```bash
    cartloader run_xenium \
      --load-xenium-ranger \
      --sge-convert \
      --run-ficture2 \
      --run-cartload2 \
      --upload-aws \
      --xenium-ranger-dir /path/to/xenium/ranger/output \
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

=== "Package SGE/cells/images"

    ```bash
    cartloader run_xenium \
      --load-xenium-ranger \
      --sge-convert \
      --import-cells \
      --import-images \
      --run-cartload2 \
      --upload-aws \
      --xenium-ranger-dir /path/to/xenium/ranger/output \
      --out-dir /path/to/out/dir \
      --s3-dir s3://example_bucket/example-id \
      --id example-id \
      --image-ids OME_DAPI \
      --spatula /path/to/spatula/binary \
      --pmtiles /path/to/pmtiles/binary \
      --tippecanoe /path/to/tippecanoe/binary \
      --aws /path/to/aws/cli/binary \
      --n-jobs 10 \
      --threads 10
    ```

<!-- === "Package existing FICTURE analysis"

    ```bash
    cartloader run_xenium \
      --load-xenium-ranger \
      --sge-convert \
      --import-cells \
      --import-images \
      --run-cartload2 \
      --upload-aws \
      --xenium-ranger-dir /path/to/xenium/ranger/output \
      --out-dir /path/to/out/dir \
      --s3-dir ${S3_PATH} \
      --id ${DATA_ID} \
      --image-ids OME_DAPI \
      --spatula /path/to/spatula/binary \
      --pmtiles /path/to/pmtiles/binary \
      --tippecanoe /path/to/tippecanoe/binary \
      --aws /path/to/aws/cli/binary \
      --n-jobs 10 \
      --threads 10
    ``` -->

---

## Parameters

### Action Parameters

See [Actions](#actions).

### Input/Output Parameters

- `--xenium-ranger-dir`: Xenium Ranger output directory (required for manual mode or with `--load-xenium-ranger`)
- `--out-dir`: Output root directory; writes `sge/`, `ficture2/`, `cartload2/` (required)
- `--xenium-ranger-assets`: Path to assets JSON from `--load-xenium-ranger` (default: `<out_dir>/xenium_ranger_assets.json`)
- `--ext-fic-dir`: External FICTURE directory to import when using `--import-ext-ficture2`

### SGE Convert (with `--sge-convert`)

- `--units-per-um`: Coordinate units per micron in raw SGE (float; default: 1.0)
- `--filter-by-density`: Enable density-based filtering (bool)
- `--exclude-feature-regex`: Regex for features to exclude (default: `^(BLANK_|DeprecatedCodeword_|NegCon|UnassignedCodeword_)`)

### FICTURE (with `--run-ficture2`)

- `--width`: Comma-separated hexagon widths in µm for train/projection (required)
- `--n-factor`: Comma-separated factor counts for training (required)
- `--colname-feature`: Feature column name (default: `gene`)
- `--colname-count`: Count column name (default: `count`)

### Import Cells (with `--import-cells`)

- `--cell-id`: Asset ID/prefix for cell outputs (default: `xeniumranger`)
- `--cell-name`: Human-readable name (defaults to `--cell-id` if omitted)
- `--tsv-cmap`: TSV colormap for cluster colors

### Import Images (with `--import-images`)

- `--image-ids`: One or more image IDs to import (default: `DAPI_OME BOUNDARY_OME INTERIOR_RNA_OME INTERIOR_PROTEIN_OME DAPI_MIP_OME`)
- `--image-colors`: HEX colors matching `--image-ids` order (defaults cover up to 10 images)
- `--all-images`: Import all images listed in `--xenium-ranger-assets`
- `--transparent-below`: Make pixels below value transparent (0–255; default: 1)

### CartLoad2 (with `--run-cartload2`)

- `--id`: Catalog ID; no whitespace; prefer `-` (required)
- `--title`: Human-readable title (quote if spaces)
- `--desc`: Short description (quote if spaces)

### Upload (with `--upload-aws` / `--upload-zenodo`)

- `--s3-dir`: S3 destination, e.g., `s3://bucket/prefix` (required for AWS)
- `--zenodo-token`: Path to file containing Zenodo token (required for Zenodo)
- `--zenodo-deposition-id`: Existing deposition ID (omit to create new)
- `--zenodo-title`: Deposition title (defaults to `--title`)
- `--creators`: One or more `"Lastname, Firstname"` entries (each quoted)

### Manual Input (enables manual mode)

- `--csv-transcript`: Transcript CSV/Parquet relative to `--xenium-ranger-dir` (used by `--sge-convert`)
- `--csv-cells`: Cell locations CSV/Parquet relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--csv-boundaries`: Cell boundary CSV relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--csv-clust`: Cluster assignments CSV relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--csv-diffexp`: Differential expression CSV relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--ome-tifs`: One or more OME‑TIFF paths relative to `--xenium-ranger-dir` (used by `--import-images`)

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

- `--load-xenium-ranger`: A Xenium assets JSON summarizing detected inputs (path: `--xenium-ranger-assets`). See example structure in the Xenium pipeline vignette.
- `--sge-convert`: Unified SGE files and an assets manifest (`sge_assets.json`). See details in [sge_convert reference](./sge_convert.md#output)
- `--run-ficture2` or `--import-ext-ficture2`: Generates FICTURE analysis artifacts (factor models, decoded maps, PMTiles, summaries). See details in [run_ficture2 reference](./run_ficture2.md#output)
- `--import-cells`: Cells PMTiles, boundaries GeoJSON, and summaries under `cartload2/`. See details in [import_cells reference](./import_cell.md#output).
- `--import-images`: Background image PMTiles (and intermediates) under `cartload2/`. See details in [import_images reference](./import_image.md#output).
- `--run-cartload2`: Sources packaged into PMTiles and a `catalog.yaml`. See details in [run_cartload2 reference](./run_cartload2.md#output).
- `--upload-aws`: Uploads generated assets and `catalog.yaml` to the specified S3 path.
- `--upload-zenodo`: Uploads `catalog.yaml` and selected assets to Zenodo.

See also example output in [Xenium end-to-end tutorial](../vignettes/pipelines/xenium.md#outputs)
