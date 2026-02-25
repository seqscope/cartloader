# SGE Add-on Modules

## Overview

`CartLoader` provides add-on modules to support specialized SGE preparation tasks that complement `sge_convert`:

- [`sge_orientate`](#reorient-sge-rotateflip): Rotate/flip SGE coordinates (e.g., to match image orientation).
- [`sge_stitch`](#stitch-multiple-sge-datasets): Stitch two or more SGE datasets into one.

---
## Requirements

- Input SGE(s) in the unified format (from [SGE format conversion](./sge_convert.md))
- Pre-installed tools: `gzip`, `spatula`, `gdal_translate`, `gdalwarp`

---
## Actions

### Reorient SGE (`sge_orientate`)
Rotations are clockwise; flips are applied after rotation.

### Stitch Multiple SGE Datasets (`sge_stitch`)
Stitch two or more SGE datasets into a single SGE, with optional per‑tile orientation. Can also run density‑based filtering and quick visualization.

---
## Example Usage

### Reorient SGE
Below is an example with 90° rotation and a horizontal flip.

```bash
OUT_PREFIX=r90_hflip                # Replace r90_hflip with your prefix
cartloader sge_orientate  \
    --in-transcript /path/to/transcripts.unsorted.tsv.gz \
    --in-feature  /path/to/feature.clean.tsv.gz \
    --in-minmax /path/to/coordinate_minmax.tsv \
    --out-transcript /path/to/${OUT_PREFIX}.transcripts.unsorted.tsv.gz \
    --out-feature /path/to/${OUT_PREFIX}.feature.clean.tsv.gz \
    --out-minmax  /path/to/${OUT_PREFIX}.coordinate_minmax.tsv \
    --rotate 90 \
    --flip-horizontal
```

### Stitch Multiple SGE Datasets
Define each input tile as a CSV of paths and tile metadata:

```bash
cartloader sge_stitch \
--in-tiles \
    "/path/to/tileA/transcripts.unsorted.tsv.gz,/path/to/tileA/feature.clean.tsv.gz,/path/to/tileA/coordinate_minmax.tsv,0,0,0,false,false" \
    "/path/to/tileB/transcripts.unsorted.tsv.gz,/path/to/tileB/feature.clean.tsv.gz,/path/to/tileB/coordinate_minmax.tsv,0,1,90,false,true" \
--out-dir /path/to/output/dir \
--units-per-um 1.0 \
--precision 2
```

---
## Parameters

### Reorient SGE

#### Input/Output Parameters
- `--in-transcript` (str): Path to input transcript SGE TSV (gzipped).
- `--in-feature` (str): Path to input feature TSV (gzipped).
- `--in-minmax` (str): Path to input coordinate min/max TSV.
- `--out-transcript` (str): Output path for oriented transcript SGE TSV (gzipped).
- `--out-feature` (str): Output path for feature TSV (copied from input).
- `--out-minmax` (str): Output path for oriented min/max TSV.

#### Key Parameters
- `--rotate` (90|180|270): Rotate clockwise (applied before flips).
- `--flip-vertical` (flag): Flip vertically (flipped along the Y‑axis; applied after rotation).
- `--flip-horizontal` (flag): Flip horizontally (flipped along the X‑axis; applied after rotation).
- `--chunk-size` (int, default: 1_000_000): Rows per chunk when processing the transcript TSV.

### Stitch Multiple SGE Datasets

#### SGE Stitch
- `--out-dir` (str): Output directory for stitched (filtered) SGE, visualizations, and Makefile.
- `--in-tiles` (list of str): One or more tile tuples describing each input tile.
    
!!! info "Input Tile Tuple Format"

    Each tile tuple should provide in the following format. Recommend to quote each tuple
    ```
    "transcript_path,feature_path,minmax_path,row,col,rotate,vertical_flip,horizontal_flip"
    ```
    
    - `transcript_path` (str): Path to the tile’s transcript‑indexed SGE TSV.gz.
    - `feature_path` (str): Path to the tile’s per‑gene UMI count TSV.gz.
    - `minmax_path` (str): Path to the tile’s coordinate min/max TSV.
    - `row`, `col` (int): Integer grid indices (0‑based) locating the tile within the stitched mosaic (row = tile grid row, col = tile grid column).
    - `rotate` (int): One of `0`, `90`, `180`, `270` (clockwise rotation).
    - `vertical_flip`, `horizontal_flip`: Boolean flags `true`/`false` (case‑insensitive); flips are applied after rotation.

??? note "Auxiliary SGE Stitching Parameters"
    We recommend using default values; override only if needed.

    **Auxiliary Output Parameters**:

    - `--out-transcript` (str): Output transcript‑indexed SGE under `--out-dir` (default: `transcripts.unsorted.tsv.gz`).
    - `--out-feature` (str): Output per‑gene UMI counts under `--out-dir` (default: `feature.clean.tsv.gz`).
    - `--out-minmax` (str): Output coordinate min/max TSV under `--out-dir` (default: `coordinate_minmax.tsv`).
    - `--out-tile-minmax` (str): Output per‑tile min/max TSV (default: `coordinate_minmax_per_tile.tsv`).
    - `--out-json` (str): Output JSON summarizing SGE assets (default: `<out-dir>/sge_assets.json`).

    **Auxiliary Colname Parameters**:

    * `--colname-count` (str): Column name of the UMI count (default: count).
    * `--colname-feature-name` (str): Column name of feature name, typically gene (default: gene).
    * `--colname-x` (str): Column name of X coordinates (default: X).
    * `--colname-y` (str): Column name of Y coordinates (default: Y).
    * `--colnames-others` (list): Column names of other information to keep.

    **Auxiliary Scale and Precision Parameters**:
    
    - `--units-per-um` (float): Coordinate units per µm for global transform (default: 1.0). If input SGE is generated by `sge_convert`, use the default value.
    - `--precision` (int): Decimal precision for stitched coordinates (default: 2).

#### Density‑based Filtering

* `--filter-by-density` (flag): Enable SGE filtering by density.
* `--out-filtered-prefix` (str): Prefix for output filtered SGE files (default: filtered).

??? note "Auxiliary Density‑based Filtering Parameters"
    We recommend using default values; override only if needed.

    * `--radius` (int): Radius for the polygon area calculation (default: 15).
    * `--quartile` (int): Quartile for the polygon area calculation (default: 2).
    * `--hex-n-move` (int): Sliding step (default: 1).
    * `--polygon-min-size` (int): Minimum polygon size (default: 500).

#### SGE Visualization
* `--sge-visual` (flag): Enable SGE visualization.
* `--north-up` (flag): Enable north‑up orientation for the SGE visualization.

??? note "Auxiliary SGE Visualization Parameters"
    We recommend using default values; override only if needed.
    
    * `--out-xy` (str): File name for output SGE visualization image (default: `xy.png`).
    * `--out-northup-tif` (str): File name for output north‑up oriented image (default: `xy_northup.tif`).
    * `--srs` (str): Spatial reference system (used with `--north-up`; default: EPSG:3857).
    * `--resample` (str): Resampling method (used with `--north-up`; options: near, bilinear, cubic, etc.; default: cubic).
    * `--gdal_translate` (str): Path to `gdal_translate` binary (used with `--north-up`; default: `gdal_translate`)
    * `--gdalwarp` (str): Path to `gdalwarp` binary (used with `--north-up`; default: `gdalwarp`)

#### Environment Parameters
* `--gzip` (str): Path to gzip binary (default: gzip).
* `--spatula` (str): Path to spatula binary (default: spatula).

---
## Output

### Reorient SGE
An output SGE matrix after orientation.

### Stitch Multiple SGE Datasets
An SGE matrix consisting of all input SGE in a specific layout. In addition:

* If `--filter-by-density` is enabled, a polygon-based filtered SGE matrix is created from the stitched SGE.
* If `--sge-visual` is enabled, a monochrome PNG image is generated to visualize the stitched SGE.
* If `--north-up` is enabled, a georeferenced TIFF image is produced with a north-up orientation.
