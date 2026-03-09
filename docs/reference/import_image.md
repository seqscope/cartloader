# Background Image Import

## Overview

`import_image` loads a histology or background image and produces a web‑ready PMTiles layer.
It supports OME‑TIFF/TIFF/PNG inputs and offers optional orientation transforms.

!!! question
    - [How do I decide orientation for background image import?](../faq/background_image_orientation.md)

---
## Requirements
- Input image: OME‑TIFF/TIFF/PNG.
- Georeference bounds (required if input image is neither georeferenced nor OME‑TIFF)
- Pre-installed CLI tools: `pmtiles`, `gdal_translate`, `gdaladdo`, `gdalinfo`.

---
## Actions
!!! warning "Action Specifications"
    No action runs by default. Activate at least one using the [action parameters](#action-parameters).

### Georeference Step (`--georeference`)
If input is OME-TIFF, extracts a PNG and bounds from the OME‑TIFF and georeferencing is applied automatically.

If `--georeference` is applied with `--georef-*`, georeferencing is applied with bounds provided via one of `--georef-*`.

### Orientation Step (`--rotate`, `--flip-vertical`, `--flip-horizontal`)
If any of `--rotate`, `--flip-vertical`, `--flip-horizontal` is set, rotation and flips can be applied accordingly prior to tiling. It supports rotate the image clockwise by 90/180/270 degrees. Rotation is applied before flips.

### OME‑TIFF to PNG Step (`--ome2png`)
Converts an OME‑TIFF to PNG, reads bounds from OME metadata, and records color mode when applicable.

### PNG to PMTiles Step (`--png2pmtiles`)
Converts the PNG to GeoTIFF and then to PMTiles. An asset JSON is written during conversion.

### Catalog Update Step (`--update-catalog`)
Appends the generated `<img-id>.pmtiles` to `catalog.yaml` as a basemap layer. This only applies if `catalog.yaml` exists.

---
## Example Usage

### 1) OME‑TIFF → PNG → PMTiles (auto georeference)

```bash
IMAGE_ID=TEST_ID      # replace TEST_ID with your own ID (no whitespace)

cartloader import_image \
  --ome2png --png2pmtiles --update-catalog \
  --in-img /path/to/input.ome.tif \
  --out-dir /path/to/out \
  --img-id  ${IMAGE_ID} \
  --pmtiles /path/to/pmtiles \
  --gdal_translate /path/to/gdal_translate \
  --gdaladdo /path/to/gdaladdo
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

# C) Pixel TSV from run_ficture2 (e.g., *.pixel.sorted.tsv.gz)
cartloader import_image \
  --png2pmtiles --georeference \
  --in-img /path/to/image.png \
  --georef-pixel-tsv /path/to/decode.pixel.sorted.tsv.gz \
  --out-dir /path/to/out \
  --img-id ${IMAGE_ID}
```

If your input image needs orientation, apply `--rotate`, `--flip-horizontal`, and/or `--flip-vertical`. For example:

```bash
cartloader import_image \
  --png2pmtiles --georeference \
  --rotate 90 --flip-horizontal \
  --in-img /path/to/image.png \
  --georef-bounds "ulx,uly,lrx,lry" \
  --out-dir /path/to/out \
  --img-id ${IMAGE_ID}
```

---
## Parameters

Below are the core arguments you’ll typically set. For all other options, expand the collapsible "Auxiliary Parameters" section.

### Action Parameters

- `--ome2png` (flag): Extract PNG (and bounds) from OME‑TIFF (auto‑georeference), and/or
- `--png2pmtiles` (flag): Convert PNG to PMTiles (via GeoTIFF/MBTiles).
- `--georeference` (flag): For manual georeference (non‑OME inputs) (used with one of `--georef-*`).
- `--rotate` (int) [90|180|270]: Rotate clockwise (applied before flips).
- `--flip-vertical` (flag): Flip vertically (around X axis; after rotation).
- `--flip-horizontal` (flag): Flip horizontally (around Y axis; after rotation).
- `--update-catalog` (flag): Update or create `catalog.yaml` with the generated PMTiles as a basemap.

### Input/Output Parameters

- `--out-dir` (str): Output directory for generated files.
- `--img-id` (str): Image ID used as filename prefix and basemap ID (when updating catalog). No whitespace.
- Input source (choose one):
    - `--in-img` (str): Path to input image (PNG or OME‑TIFF/TIFF).
    - `--in-json` (str): JSON/YAML mapping from image IDs to file paths.

??? note "Auxiliary Parameters"

    **Auxiliary for --ome2png**

      - `--micron2pixel-csv` (str): Vizgen transform CSV (micron_to_mosaic_pixel_transform.csv).
      - `--page` (int): Z‑slice index to extract from multi‑page OME‑TIFF (3D).
      - `--level` (int): Resolution level index to extract from OME‑TIFF.
      - `--series` (int): Series index to extract from OME‑TIFF.
      - `--upper-thres-quantile` / `--upper-thres-intensity`: Rescale cap (mutually exclusive).
      - `--lower-thres-quantile` / `--lower-thres-intensity`: Rescale floor (mutually exclusive).
      - `--transparent-below` (int): Make pixels below threshold transparent (0–255).
      - `--colorize` (str): Colorize mono images using an RGB hex or name.
      - `--high-memory` (flag): Use a memory‑intensive path for large OME‑TIFFs.

    **Auxiliary for --png2pmtiles**

      - `--srs` (str, default: `EPSG:3857`): Spatial reference.
      - `--mono` (flag): Input PNG is single‑band mono (skip if `--ome2png`).
      - `--rgba` (flag): Input PNG is 4‑band RGBA (skip if `--ome2png`).
      - `--resample` (str, default: `cubic`): GDAL resampling.
      - `--blocksize` (int, default: `512`): GDAL block size in pixels.

    **Auxiliary for --georeference**
      Pick one bounds source:
      - `--georef-detect` (flag): Used the detect bounds from image metadata (e.g., 'OME'). Extracted bounds will automatically be applied to png2pmtiles.
      - `--georef-pixel-tsv` (str): Pixel TSV (e.g., `*.pixel.sorted.tsv.gz` from run_ficture2).
      - `--georef-bounds-tsv` (str): One‑line TSV with `ulx,uly,lrx,lry`.
      - `--georef-bounds` (str): Bounds string `"ulx,uly,lrx,lry"`.

    **Environment Parameters**

      - `--pmtiles` (str, default: `pmtiles`): Path to go‑pmtiles binary.
      - `--gdal_translate` (str, default: `gdal_translate`)
      - `--gdaladdo` (str, default: `gdaladdo`)
      - `--gdalinfo` (str, default: `gdalinfo`)

    **Run Parameters**

      - `--dry-run` (flag): Generate the Makefile; do not execute.
      - `--restart` (flag): Ignore existing outputs and rerun steps.
      - `--n-jobs` (int): Number of parallel jobs (default: 1).
      - `--makefn` (str): Makefile name to write (default: `<out-dir>/<img-id>.mk`).

---
## Output

Outputs are written under `--out-dir` with prefix `<img-id>`:

- `<img-id>.pmtiles`: Final PMTiles for web visualization.
- `<img-id>.json`: Asset JSON produced during PNG→PMTiles conversion.
- `<img-id>.geotif.tif` and intermediate tiles (when relevant).
