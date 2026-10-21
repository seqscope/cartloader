# Visium HD Pipeline Orchestrator (Legacy)

!!! warning "Legacy"
    `run_visiumhd` is superseded by [`run_together`](./run_together.md) (see the [Visium HD profile](./platforms/visium_hd.md)). It is kept for reference and is not part of the recommended workflow; see [Legacy Orchestrators](./legacy.md) for what `run_together` covers and the remaining gaps.

    **Known issues:** combining `--run-ficture2` with `--run-cartload2` fails, because the packaging step waits for decode outputs named `t<w>_f<k>_p<w>_a6.*` that `run_ficture2` no longer writes (they are now `t<w>_f<k>_pixel.*`); `--import-ext-ficture2` fails, because it passes `--ext-fic-dir`, which `run_cartload2` does not accept; and `--run-cartload2` without FICTURE fails as well. Run [`run_cartload2`](./run_cartload2.md) directly instead. See the [full list and workarounds](./legacy.md#known-issues).

## Overview

An end‑to‑end workflow for one 10x Visium HD sample that runs `CartLoader` steps (detect Space Ranger assets, convert the SGE, run FICTURE, import cells, square bins, and background images, package assets, and upload) as needed.

---
## Requirements

- Space Ranger output
- Tools depending on the actions: `spatula`, `punkst`, `tippecanoe`, `pmtiles`, `geotiff2pmtiles` (bundled submodules by default), plus `gzip`, `sort`, `python3`, GDAL, `pigz` or `parquet-tools`, `Rscript`, and the `aws` CLI on `PATH`

---
## Actions

!!! info "Action Specifications"

    `run_visiumhd` orchestrates multiple `CartLoader` modules; enable and combine actions with flags.

    No action runs by default. Activate at least one action flag.

See details for each flag:

- `--load-space-ranger`: Detect Space Ranger outputs and write the assets JSON ([`load_space_ranger`](./visiumhd_addon.md)); archives are extracted to `<out-dir>/unzip`. Skipped if the assets JSON exists, unless `--restart`.
- [`--sge-convert`](./sge_convert.md#actions)
- [`--run-ficture2`](./run_ficture2.md#actions): Runs with `--decode-scale 2`.
- `--import-ext-ficture2`: Package an existing FICTURE directory (`--ext-fic-dir`); currently fails (see above).
- [`--import-cells`](./import_cell.md#actions): Runs `import_visiumhd_cell --all`; skipped if `<out-dir>/cartload2/<cell-id>_assets.json` exists, unless `--restart`.
- [`--import-squares`](./import_square.md#actions): Runs `import_visiumhd_square` per bin size, writing `<square-id>-sq<NNN>_assets.json` (e.g., `spaceranger-sq008_assets.json`). In auto-detect mode, every `GRID_<N>um` manifest entry with `POSITION`, `CLUSTER`, `DE`, and `UMAP_PROJ` is imported; incomplete entries are skipped.
- [`--import-images`](./import_image.md#actions): Imports each image with `--georef-detect ome`, so the image must carry OME-XML pixel-size metadata (e.g., a Space Ranger H&E BTF).
- [`--run-cartload2`](./run_cartload2.md#actions)
- [`--upload-aws`](./upload_aws.md#actions)
- [`--upload-zenodo`](./upload_zenodo.md#actions): Runs `upload_zenodo --upload-method catalog` (with `--overwrite`, not `--restart`, when `--restart` is set); skipped when its `.zenodo.done` flags exist.

`--import-cells`, `--import-squares`, and `--import-images` also update `catalog.yaml` when it already exists and `--run-cartload2` is not set.

---
## Example Usage

Tool paths default to the bundled submodules; add `--spatula`, `--ficture2`, `--pmtiles`, or `--tippecanoe` only to override them.

=== "Run all modules together"

    ```bash
    cartloader run_visiumhd \
      --load-space-ranger \
      --sge-convert \
      --run-ficture2 \
      --import-cells \
      --import-squares \
      --import-images \
      --run-cartload2 \
      --upload-aws \
      --space-ranger-dir /path/to/space/ranger/output \
      --out-dir /path/to/out/dir \
      --s3-dir s3://example_bucket/example-id \
      --width 12 \
      --n-factor 12 \
      --id example-id \
      --image-ids HnE \
      --n-jobs 10 \
      --threads 10
    ```

---
## Parameters

Below are the core arguments you’ll typically set. For all other options, expand the collapsible "Auxiliary Parameters" section.

### Action Parameters

See [Actions](#actions).

### Input/Output Parameters

- `--space-ranger-dir`: Space Ranger output directory (required for manual mode or with `--load-space-ranger`)
- `--out-dir`: Output root directory; writes `sge/`, `ficture2/`, `cartload2/` (required)
- `--space-ranger-assets`: Path to the assets JSON from `--load-space-ranger` (default: `<out_dir>/space_ranger_assets.json`)
- `--ext-fic-dir`: External FICTURE directory to import with `--import-ext-ficture2`

### SGE Convert (with `--sge-convert`)

- `--units-per-um`: Coordinate units per micron in the raw SGE (float; manual mode only, when no `--json-scale` is given)
- `--exclude-feature-regex`: Regex for features to exclude
- `--filter-by-density`: Enable density-based filtering (bool)

In auto-detect mode (with `--load-space-ranger` or an existing assets JSON), SGE inputs are read from `--space-ranger-assets`.
In manual mode, provide inputs relative to `--space-ranger-dir`:

- `--mex-transcript`: Transcript MEX directory (e.g., `binned_outputs/square_002um/filtered_feature_bc_matrix`)
- `--parquet-position`: Positions Parquet (e.g., `binned_outputs/square_002um/spatial/tissue_positions.parquet`)
- `--json-scale` or `--units-per-um`: Scale source (one is required)

### FICTURE (with `--run-ficture2`)
- `--width`: Comma‑separated hexagon widths in µm for train/projection (required)
- `--n-factor`: Comma‑separated factor counts for training (required)
- `--colname-feature`: Feature column name (default: `gene`)
- `--colname-count`: Count column name (default: `count`)
- `--fic-include-feature-regex` (str): Regex of feature names to include in FICTURE analysis.
- `--fic-exclude-feature-regex` (str): Regex of feature names to exclude in FICTURE analysis.

### Import Cells (with `--import-cells`)

- `--cell-id`: Asset ID/prefix for cell outputs (default: `spaceranger`)
- `--cell-name`: Human‑readable name (defaults to `--cell-id`)
- `--tsv-cmap`: TSV colormap for cluster colors (also used by `--import-squares`)

Manual mode inputs (relative to `--space-ranger-dir`):

- `--geojson-cells`: Cell segmentation GeoJSON
- `--mtx-cells`: Cell-feature MEX directory
- `--csv-clust`: Cluster assignments CSV
- `--csv-diffexp`: Differential expression CSV

### Import Squares (with `--import-squares`)

- `--square-id`: Prefix for square-bin assets (default: `spaceranger`)
- `--use-parquet-tools`: Use `parquet-tools` instead of polars/pigz for Parquet→CSV conversion (slower on large files)
- `--square-input`: Manual mode only: one or more square bin definitions, each as `bin_size,bin_pos_parquet,bin_csv_cluster,bin_csv_diffexp,bin_scale_json,bin_csv_umap` (paths relative to `--space-ranger-dir`). The UMAP CSV is required by `import_visiumhd_square`; see [known issues](./legacy.md#known-issues) for the bin MEX.

### Import Images (with `--import-images`)

- `--image-ids`: One or more image IDs to import (default: `HnE`)
- `--all-images`: Import all images listed in `--space-ranger-assets`

Manual mode inputs (relative to `--space-ranger-dir`):

- `--tifs`: One or more BTF/TIFF paths, in `--image-ids` order

### CartLoad2 (with `--run-cartload2`)

- `--id`: Catalog ID; no whitespace; prefer `-` over `_` (required)
- `--title`: Human‑readable title (quote if spaces)
- `--desc`: Short description (quote if spaces)

### Upload (with `--upload-aws` / `--upload-zenodo`)

- `--s3-dir`: S3 destination, e.g., `s3://bucket/prefix` (required for AWS)
- `--profile`: AWS CLI profile for the upload
- `--zenodo-token`: Path to a file containing the Zenodo token (required for Zenodo)
- `--zenodo-deposition-id`: Existing deposition ID (omit to create new)
- `--zenodo-title`: Deposition title (defaults to `--title`)
- `--creators`: One or more `"Lastname, Firstname"` entries (each quoted; required for a new deposition)

??? info "Auxiliary Parameters"

    #### Run Options

    - `--dry-run`: Print commands without executing (bool; default: false)
    - `--restart`: Ignore existing outputs and rerun steps (bool; default: false)
    - `--n-jobs`: Parallel jobs to run (int; default: 1)
    - `--threads`: Max threads per job, e.g., tippecanoe (int; default: 4)

    #### Environment / Tools

    - `--pmtiles`: Path to `pmtiles` (default: `<cartloader_dir>/submodules/pmtiles/pmtiles`)
    - `--tippecanoe`: Path to `tippecanoe` (default: `<cartloader_dir>/submodules/tippecanoe/tippecanoe`)
    - `--spatula`: Path to `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`)
    - `--ficture2`: Path to the punkst repository (default: `<cartloader_dir>/submodules/punkst`)
    - `--gdal_translate`, `--gdaladdo`, `--gdalwarp`: GDAL binaries (default: on `PATH`)
    - `--parquet-tools`: Path to `parquet-tools` (default: each command's own default, `parquet-tools`)
    - `--pigz`, `--pigz-threads`: pigz binary and threads for Parquet conversion (default: `pigz`, `4`)
    - `--aws`: Path to the AWS CLI (default: `aws`)
    - `--R`, `--gzip`, `--sort`, `--python`: Binaries on `PATH` (default: `Rscript`, `gzip`, `sort`, `python3`)

---
## Output

Outputs depend on the enabled action flags. Refer to each module page for full output details and file formats.

- `--load-space-ranger`: A Space Ranger assets JSON summarizing detected inputs (path: `--space-ranger-assets`). See the example in [Visium HD Add-ons](./visiumhd_addon.md#output).
- `--sge-convert`: Unified SGE files and an assets manifest (`sge_assets.json`) under `sge/`. See [sge_convert reference](./sge_convert.md#output).
- `--run-ficture2`: FICTURE analysis artifacts under `ficture2/`. See [run_ficture2 reference](./run_ficture2.md#output).
- `--import-cells`: Cell and boundary PMTiles and a cell assets JSON under `cartload2/`. See [import_cell reference](./import_cell.md#output).
- `--import-squares`: Square-bin PMTiles plus one assets JSON per bin size (e.g., `<square-id>-sq008_assets.json`) under `cartload2/`. See [import_square reference](./import_square.md#output).
- `--import-images`: Background image PMTiles (and intermediates) under `cartload2/`. See [import_image reference](./import_image.md#output).
- `--run-cartload2`: PMTiles and a `catalog.yaml` under `cartload2/`. See [run_cartload2 reference](./run_cartload2.md#output).
- `--upload-aws`: Uploads the catalog-listed assets and `catalog.yaml` to `--s3-dir`.
- `--upload-zenodo`: Uploads `catalog.yaml` and catalog-listed assets to Zenodo.
