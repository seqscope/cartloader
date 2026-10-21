# Background Image Import

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) calls `import_image` in its `images` stage for each image modality you declare (see [Image Modalities](./run_together_images.md)). This page is for importing an image directly, step by step.

## Overview

`import_image` loads a histology or background image (OME‑TIFF, TIFF, or PNG) and produces a georeferenced, web‑ready raster PMTiles layer, with optional orientation transforms and catalog registration.

!!! question
    - [How do I decide orientation for background image import?](../faq/background_image_orientation.md)

---
## Requirements
- Input image: OME‑TIFF, TIFF, or PNG.
- A bounds source unless the image already carries a geotransform: OME metadata (`--ome2png` or `--georef-detect ome`), or explicit bounds (`--georef-bounds`, `--georef-bounds-tsv`, `--georef-pixel-tsv`).
- Tools: `geotiff2pmtiles` (default: bundled submodule), and GDAL's `gdalwarp`, `gdalinfo`, and `gdal_translate` on `PATH` (`gdal_translate` is needed for georeferencing).
- For `--update-catalog`: an existing `catalog.yaml` (e.g., from [`run_cartload2`](./run_cartload2.md)).

---
## Actions
!!! warning "Action Specifications"
    No action runs by default. Pass at least one of `--ome2png`, `--png2pmtiles`, `--georeference`, `--rotate`, `--flip-vertical`, or `--flip-horizontal`. Georeferencing and orientation are applied inside the PNG‑to‑PMTiles step, so they take effect only together with `--png2pmtiles`.

### OME‑TIFF to PNG Step (`--ome2png`)
Extracts a PNG from an OME‑TIFF (`cartloader image_ome2png`), reads its bounds from the OME metadata (or from `--micron2pixel-csv` / `--px-per-um-x/y`), and records its color mode. This also turns on `--georeference` with those bounds, so `--georef-bounds`, `--georef-bounds-tsv`, and `--georef-pixel-tsv` cannot be combined with it.

### Georeference Step (`--georeference`)
Assigns bounds and the spatial reference (`--srs`) to the image with `gdal_translate`. Provide exactly one bounds source: `--georef-bounds`, `--georef-bounds-tsv`, `--georef-pixel-tsv`, or `--georef-detect` (`ome` or `gtiff`).

### Orientation Step (`--rotate`, `--flip-vertical`, `--flip-horizontal`)
Rotates the georeferenced image clockwise by 90/180/270 degrees and/or flips it (rotation before flips) with `gdalwarp`. Equivalent combinations are normalized first (e.g., `--rotate 180` becomes both flips).

### PNG to PMTiles Step (`--png2pmtiles`)
Converts the (georeferenced, oriented) image into `<img-id>.pmtiles` with `geotiff2pmtiles` and writes `<img-id>_assets.json`. The image must carry a coordinate system at this point; otherwise the step stops with an error rather than writing an empty layer.

### Catalog Update Step (`--update-catalog`)
Adds `<img-id>.pmtiles` to the catalog as a basemap with ID `<img-id>`. The catalog (`--catalog-yaml`, default `<out-dir>/catalog.yaml`) must already exist; otherwise the run stops because make finds no rule for it.

---
## Example Usage

### 1) OME‑TIFF → PNG → PMTiles (auto georeference)

```bash
IMAGE_ID=TEST_ID      # replace TEST_ID with your own ID (no whitespace)

cartloader import_image \
  --ome2png --png2pmtiles --update-catalog \
  --in-img /path/to/input.ome.tif \
  --out-dir /path/to/out \
  --img-id  ${IMAGE_ID}
```

### 2) PNG → PMTiles (manual georeference via bounds)

Pick one bounds input method and enable `--georeference`.

```bash
IMAGE_ID=TEST_ID      # replace TEST_ID with your own ID (no whitespace)

# A) Bounds string
cartloader import_image \
  --png2pmtiles --georeference \
  --in-img /path/to/image.png \
  --georef-bounds "ulx,uly,lrx,lry" \
  --out-dir /path/to/out \
  --img-id ${IMAGE_ID}

# B) Bounds TSV
cartloader import_image \
  --png2pmtiles --georeference \
  --in-img /path/to/image.png \
  --georef-bounds-tsv /path/to/bounds.tsv \
  --out-dir /path/to/out \
  --img-id ${IMAGE_ID}

# C) OME-TIFF metadata (e.g., an 8-bit RGB H&E OME-TIFF), without extracting a PNG first
cartloader import_image \
  --png2pmtiles --georeference \
  --in-img /path/to/image.ome.tif \
  --georef-detect ome \
  --out-dir /path/to/out \
  --img-id ${IMAGE_ID}
```

If your input image needs orientation, add `--rotate`, `--flip-horizontal`, and/or `--flip-vertical`. For example:

```bash
cartloader import_image \
  --png2pmtiles --georeference \
  --rotate 90 --flip-horizontal \
  --in-img /path/to/image.png \
  --georef-bounds "ulx,uly,lrx,lry" \
  --out-dir /path/to/out \
  --img-id ${IMAGE_ID}
```

!!! warning "`--rotate 90/270` with `--ome2png`"
    Combining `--ome2png` with a rotation that normalizes to 90 degrees (`--rotate 90` or `--rotate 270`, alone or with flips) currently fails. Extract the PNG first, or rotate in a separate `--png2pmtiles` run on the extracted PNG. `--rotate 180` and flips work with `--ome2png`.

---
## Parameters

Below are the core arguments you’ll typically set. For all other options, expand the collapsible "Auxiliary Parameters" section.

### Action Parameters

- `--ome2png` (flag): Extract a PNG and its bounds from an OME‑TIFF (enables georeferencing).
- `--png2pmtiles` (flag): Convert the image to PMTiles.
- `--georeference` (flag): Georeference the image using one `--georef-*` source.
- `--rotate` (str) [90|180|270]: Rotate clockwise (applied before flips).
- `--flip-vertical` (flag): Flip vertically (around the X axis; after rotation).
- `--flip-horizontal` (flag): Flip horizontally (around the Y axis; after rotation).
- `--update-catalog` (flag): Add the PMTiles to an existing catalog as a basemap.
- `--skip-image-errors` (flag): With `--ome2png`, skip an unreadable or corrupt OME‑TIFF with a warning (no PMTiles, no catalog entry) instead of failing.

### Input/Output Parameters

- `--out-dir` (str, required): Output directory for generated files.
- `--img-id` (str, required): Image ID used as the filename prefix and basemap ID. No whitespace.
- Input source (choose one):
    - `--in-img` (str): Path to the input image (PNG or OME‑TIFF/TIFF).
    - `--in-json` (str): JSON/YAML mapping image IDs to paths (its `IMAGES` section if present, e.g., a `load_xenium_ranger` / `load_space_ranger` manifest); `--img-id` selects the entry and `--in-img` is ignored.
- `--catalog-yaml` (str): Catalog to update with `--update-catalog` (default: `<out-dir>/catalog.yaml`).

??? note "Auxiliary Parameters"

    **Auxiliary for --ome2png**

      - `--micron2pixel-csv` (str): Vizgen transform CSV (`micron_to_mosaic_pixel_transform.csv`).
      - `--px-per-um-x` / `--px-per-um-y` (float): Pixel size for a plain TIFF without pixel-size metadata (both required together).
      - `--offset-px-x` / `--offset-px-y` (float): Pixel offsets used with `--px-per-um-x/y` (default: 0).
      - `--page` (int): Z‑slice index to extract from a multi‑page OME‑TIFF.
      - `--use-middle-page` (flag): Pick the middle page automatically when `--page` is not given.
      - `--level` (int): Resolution level to extract (image_ome2png default: 0).
      - `--series` (int): Series index to extract (image_ome2png default: 0).
      - `--upper-thres-quantile` / `--upper-thres-intensity`: Rescale cap (mutually exclusive).
      - `--lower-thres-quantile` / `--lower-thres-intensity`: Rescale floor (mutually exclusive).
      - `--transparent-below` (int): Make pixels below this value transparent (0–255).
      - `--colorize` (str): Colorize a mono image with an RGB hex code or color name.
      - `--high-memory` (flag): Use a faster, memory‑intensive path for large OME‑TIFFs.
      - `--shrink-factor` (float): Downsample by this factor before processing.

    **Auxiliary for --georeference**

      Pick one bounds source:

      - `--georef-bounds` (str): Bounds string `"ulx,uly,lrx,lry"`.
      - `--georef-bounds-tsv` (str): File with one line `ulx,uly,lrx,lry`.
      - `--georef-pixel-tsv` (str): A FICTURE `*.pixel.sorted.tsv.gz` whose header records the offsets and size.
      - `--georef-detect` (str): Read bounds from the image itself: `ome` (OME‑XML pixel size and offsets) or `gtiff` (keep the existing geotransform and only assign `--srs`).

      Shared:

      - `--georef-offset-x` / `--georef-offset-y` (float): Offsets added to the bounds (default: 0.0).

    **Auxiliary for --png2pmtiles**

      - `--srs` (str): Spatial reference system (default: `EPSG:3857`).
      - `--mono` / `--rgba` (flag): Declare a single-band or 4-band input (ignored with `--ome2png`, which records the color mode).
      - `--min-zoom` (int): Minimum zoom level (default: `6`).
      - `--max-zoom` (int): Maximum zoom level (default: determined by `geotiff2pmtiles`).
      - `--tile-format` (str): `png` or `webp` (default: `png`).
      - `--rescale` (str): `auto`, `linear`, `log`, or `none`; needed for 16-bit images, together with `--rescale-range "min,max"`.
      - `--resample` (str), `--blocksize` (int): Used only by the GDAL/MBTiles conversion path, which `import_image` does not select (defaults: `cubic`, `512`).
      - `--gdal-only` (flag): Currently not usable; it is forwarded to `image_png2pmtiles`, which does not accept it.

    **Environment Parameters**

      - `--geotiff2pmtiles` (str): `geotiff2pmtiles` binary (default: `<cartloader_dir>/submodules/geotiff2pmtiles/geotiff2pmtiles`).
      - `--gdal_translate` (str): `gdal_translate` binary (default: `gdal_translate`).
      - `--pmtiles` (str), `--gdaladdo` (str): Forwarded, but used only by the GDAL/MBTiles path (defaults: `<cartloader_dir>/submodules/pmtiles/pmtiles`, `gdaladdo`).
      - `--gdalinfo` (str): Accepted but not forwarded; `gdalinfo` and `gdalwarp` are taken from `PATH`.

    **Run Parameters**

      - `--dry-run` (flag): Generate the Makefile; do not execute.
      - `--restart` (flag): Ignore existing outputs and rerun steps.
      - `--n-jobs` (int): Number of parallel jobs (default: 1).
      - `--makefn` (str): Makefile name under `--out-dir` (default: `<img-id>.mk`).

---
## Output

Outputs are written under `--out-dir` with prefix `<img-id>`:

- `<img-id>.pmtiles`: Final raster PMTiles.
- `<img-id>_assets.json`: Asset JSON (`{"<img-id>": "<path to PMTiles>"}`) for `run_cartload2 --background-assets`.
- Intermediates (kept): `<img-id>.georef.tif` (with georeferencing), then either `<img-id>.warped.tif` or, with orientation, `<img-id>.<rot90|vflip|hflip...>.tif` plus a `.dim.tsv` size file.
- With `--ome2png`: `<img-id>.transform.png`, `<img-id>.transform.bounds.csv`, `<img-id>.transform.color.csv`, and a `.transform.done` flag (`.transform.skipped` when an image is skipped under `--skip-image-errors`).
- Makefiles `<img-id>.mk` and `<img-id>.png2pmtiles.mk`; `<img-id>.pmtiles.yaml.done` after a catalog update.
