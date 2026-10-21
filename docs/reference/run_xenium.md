# Xenium Pipeline Orchestrator (Legacy)

!!! warning "Legacy"
    `run_xenium` is superseded by [`run_together`](./run_together.md) (see the [Xenium profile](./platforms/xenium.md)). It is kept for reference and is not part of the recommended workflow; see [Legacy Orchestrators](./legacy.md) for what `run_together` covers and the remaining gaps.

    **Known issues:** combining `--run-ficture2` with `--run-cartload2` fails, because the packaging step waits for decode outputs named `t<w>_f<k>_p<w>_a6.*` that `run_ficture2` no longer writes (they are now `t<w>_f<k>_pixel.*`); `--import-ext-ficture2` fails, because it passes `--ext-fic-dir`, which `run_cartload2` does not accept; and `--run-cartload2` without FICTURE fails as well. Run [`run_cartload2`](./run_cartload2.md) directly instead. See the [full list and workarounds](./legacy.md#known-issues).

## Overview

An end‑to‑end workflow for one 10x Xenium sample that runs `CartLoader` steps (detect inputs, convert the SGE, run FICTURE, import cell analysis results and background images, package assets, and upload) as needed.

---
## Requirements

- Xenium Ranger output
- Tools depending on the actions: `spatula`, `punkst`, `tippecanoe`, `pmtiles`, `geotiff2pmtiles` (bundled submodules by default), plus `gzip`, `sort`, `python3`, GDAL, `parquet-tools`, `Rscript`, and the `aws` CLI on `PATH`

---
## Actions

!!! warning "Action Specifications"
    `run_xenium` orchestrates multiple `CartLoader` modules; enable and combine actions with flags.

    No action runs by default. Activate at least one action flag.

Click to see action details per flag:

* `--load-xenium-ranger`: Detect Xenium Ranger outputs and write the assets JSON ([`load_xenium_ranger`](./xenium_addon.md)); archives are extracted to `<out-dir>/unzip`. Skipped if the assets JSON exists, unless `--restart`.
* [`--sge-convert`](./sge_convert.md#actions)
* [`--run-ficture2`](./run_ficture2.md#actions)
* `--import-ext-ficture2`: Package an existing FICTURE directory (`--ext-fic-dir`); currently fails (see above).
* [`--import-cells`](./import_cell.md#actions): Runs `import_xenium_cell --all`; skipped if `<out-dir>/cartload2/<cell-id>_assets.json` exists, unless `--restart`.
* [`--import-images`](./import_image.md#actions)
* [`--run-cartload2`](./run_cartload2.md#actions)
* [`--upload-aws`](./upload_aws.md#actions)
* [`--upload-zenodo`](./upload_zenodo.md#actions): Runs `upload_zenodo --upload-method catalog` (with `--overwrite`, not `--restart`, when `--restart` is set); skipped when its `.zenodo.done` flags exist.

`--import-cells` and `--import-images` also update `catalog.yaml` when it already exists and `--run-cartload2` is not set.

---
## Example Usage

!!! warning "Replace placeholders"

    Replace placeholder paths (e.g., `/path/to/xenium/ranger/output`) or example values (e.g., `example-id` and `12`) before execution. Tool paths default to the bundled submodules; add `--spatula`, `--ficture2`, `--pmtiles`, or `--tippecanoe` only to override them.

=== "All modules (end‑to‑end)"

    The example below imports only the `DAPI_OME` image for illustration.

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
      --image-ids DAPI_OME \
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
      --image-ids DAPI_OME \
      --n-jobs 10 \
      --threads 10
    ```

---

## Parameters

### Action Parameters

See [Actions](#actions).

### Input/Output Parameters

- `--xenium-ranger-dir`: Xenium Ranger output directory (required for manual mode or with `--load-xenium-ranger`)
- `--out-dir`: Output root directory; writes `sge/`, `ficture2/`, `cartload2/` (required)
- `--xenium-ranger-assets`: Path to the assets JSON from `--load-xenium-ranger` (default: `<out_dir>/xenium_ranger_assets.json`)
- `--ext-fic-dir`: External FICTURE directory to import with `--import-ext-ficture2`

### SGE Convert (with `--sge-convert`)

- `--units-per-um`: Coordinate units per micron in the raw SGE (float; default: the `sge_convert` default, 1.0)
- `--filter-by-density`: Enable density-based filtering (bool)
- `--exclude-feature-regex`: Regex for features to exclude (default: the `sge_convert` Xenium default, `^(BLANK_|Blank-|DeprecatedCodeword_|NegCon|UnassignedCodeword_)`)

### FICTURE (with `--run-ficture2`)

- `--width`: Comma-separated hexagon widths in µm for train/projection (required)
- `--n-factor`: Comma-separated factor counts for training (required)
- `--colname-feature`: Feature column name (default: `gene`)
- `--colname-count`: Count column name (default: `count`)
- `--fic-include-feature-regex` / `--fic-exclude-feature-regex`: Regex of features to include/exclude in FICTURE

### Import Cells (with `--import-cells`)

- `--cell-id`: Asset ID/prefix for cell outputs (default: `xeniumranger`)
- `--cell-name`: Human-readable name (defaults to `--cell-id` if omitted)
- `--tsv-cmap`: TSV colormap for cluster colors

### Import Images (with `--import-images`)

- `--image-ids`: One or more image IDs to import (default: `DAPI_OME BOUNDARY_OME INTERIOR_RNA_OME INTERIOR_PROTEIN_OME DAPI_MIP_OME`; every listed ID must exist in the assets JSON)
- `--image-colors`: HEX colors matching `--image-ids` order (defaults cover up to 10 images)
- `--all-images`: Import all images listed in `--xenium-ranger-assets`
- `--transparent-below`: Make pixels below value transparent (0–255; default: 1)

### CartLoad2 (with `--run-cartload2`)

- `--id`: Catalog ID; no whitespace; prefer `-` (required)
- `--title`: Human-readable title (quote if spaces)
- `--desc`: Short description (quote if spaces)

### Upload (with `--upload-aws` / `--upload-zenodo`)

- `--s3-dir`: S3 destination, e.g., `s3://bucket/prefix` (required for AWS)
- `--profile`: AWS CLI profile for the upload
- `--zenodo-token`: Path to a file containing the Zenodo token (required for Zenodo)
- `--zenodo-deposition-id`: Existing deposition ID (omit to create new)
- `--zenodo-title`: Deposition title (defaults to `--title`)
- `--creators`: One or more `"Lastname, Firstname"` entries (each quoted; required for a new deposition)

### Manual Input (enables manual mode)

- `--csv-transcript`: Transcript CSV/Parquet relative to `--xenium-ranger-dir` (used by `--sge-convert`)
- `--csv-cells`: Cell locations CSV/Parquet relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--csv-boundaries`: Cell boundary CSV relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--csv-clust`: Cluster assignments CSV relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--csv-diffexp`: Differential expression CSV relative to `--xenium-ranger-dir` (used by `--import-cells`)
- `--ome-tifs`: One or more OME‑TIFF paths relative to `--xenium-ranger-dir`, in `--image-ids` order (used by `--import-images`)

??? info "Auxiliary Parameters"

    #### Run Options

    - `--dry-run`: Print commands without executing (bool; default: false)
    - `--restart`: Ignore existing outputs and rerun steps (bool; default: false)
    - `--n-jobs`, `-j`: Parallel jobs to run (int; default: 1)
    - `--threads`: Max threads per job, e.g., tippecanoe (int; default: 4)

    #### Environment / Tools

    - `--pmtiles`: Path to `pmtiles` (default: `<cartloader_dir>/submodules/pmtiles/pmtiles`)
    - `--tippecanoe`: Path to `tippecanoe` (default: `<cartloader_dir>/submodules/tippecanoe/tippecanoe`)
    - `--spatula`: Path to `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`)
    - `--ficture2`: Path to the punkst repository (default: `<cartloader_dir>/submodules/punkst`)
    - `--gdal_translate`, `--gdaladdo`, `--gdalwarp`: GDAL binaries (default: on `PATH`)
    - `--parquet-tools`: Path to `parquet-tools` (default: each command's own default, `parquet-tools`)
    - `--aws`: Path to the AWS CLI (default: `aws`)
    - `--R`, `--gzip`, `--sort`, `--python`: Binaries on `PATH` (default: `Rscript`, `gzip`, `sort`, `python3`)

---
## Output

Outputs depend on the enabled action flags. Refer to each module page for full output details and file formats.

- `--load-xenium-ranger`: A Xenium assets JSON summarizing detected inputs (path: `--xenium-ranger-assets`). See the example in [Xenium Add-ons](./xenium_addon.md#output).
- `--sge-convert`: Unified SGE files and an assets manifest (`sge_assets.json`) under `sge/`. See [sge_convert reference](./sge_convert.md#output).
- `--run-ficture2`: FICTURE analysis artifacts under `ficture2/`. See [run_ficture2 reference](./run_ficture2.md#output).
- `--import-cells`: Cell and boundary PMTiles and a cell assets JSON under `cartload2/`. See [import_cell reference](./import_cell.md#output).
- `--import-images`: Background image PMTiles (and intermediates) under `cartload2/`. See [import_image reference](./import_image.md#output).
- `--run-cartload2`: PMTiles and a `catalog.yaml` under `cartload2/`. See [run_cartload2 reference](./run_cartload2.md#output).
- `--upload-aws`: Uploads the catalog-listed assets and `catalog.yaml` to `--s3-dir`.
- `--upload-zenodo`: Uploads `catalog.yaml` and catalog-listed assets to Zenodo.

See also the example output in the [Xenium end-to-end tutorial](../vignettes/pipelines/xenium.md#outputs).
