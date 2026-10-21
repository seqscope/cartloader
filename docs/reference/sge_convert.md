# Spatial Gene Expression Format Conversion

## Overview

SGE datasets vary widely in format and resolution across platforms. Since `FICTURE` requires SGE in a specific format, `CartLoader` provides the `sge_convert` module to standardize raw, platform-specific SGE into a FICTURE-compatible format (coordinates in µm).

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) calls `sge_convert` as its ingest stage, with per-platform settings. Use this page to run it directly in step-by-step or custom workflows.

!!! question
    - [How do I define scaling (`--units-per-um`) for `sge_convert`?](../faq/define_units_per_um.md)

---
## Requirements

The input (raw, platform‑specific SGE) must provide at least:

* Spatial coordinates (X, Y)
* Feature names (such as gene symbols)
* Expression counts (optional for CSV/TSV input; each row counts as 1 when no count column is given)

!!! info "Platform Compatibility"

    | Source                                                                            | `--platform`      | Input flags                                                                 |
    |-----------------------------------------------------------------------------------|-------------------|-----------------------------------------------------------------------------|
    | [10x Visium HD](https://www.10xgenomics.com/platforms/visium)                     | `10x_visium_hd`   | `--in-mex` and `--pos-parquet` (required); `--scale-json` (recommended)     |
    | [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0)                   | `seqscope`        | `--in-mex`, **or** `--in-csv` (raw per-molecule TSV)                        |
    | Illumina StrataMap                                                                | `illumina`        | `--in-mex` (coordinates are embedded in the barcodes)                       |
    | [10x Xenium](https://www.10xgenomics.com/platforms/xenium)                        | `10x_xenium`      | `--in-csv` **or** `--in-parquet` (not both)                                 |
    | [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq) | `bgi_stereoseq`   | `--in-csv`                                                                  |
    | [CosMx SMI](https://nanostring.com/products/cosmx-spatial-molecular-imager)       | `cosmx_smi`       | `--in-csv`                                                                  |
    | [Vizgen MERSCOPE](https://vizgen.com/merscope-ultra/)                             | `vizgen_merscope` | `--in-csv`                                                                  |
    | [Pixel-seq](https://www.cell.com/cell/fulltext/S0092-8674(22)01367-8)             | `pixel_seq`       | `--in-csv`                                                                  |
    | Nova-ST                                                                           | `nova_st`         | `--in-csv`                                                                  |
    | Generic CSV/TSV input <sup>1</sup>                                                | `generic`         | `--in-csv`                                                                  |

    For `10x_xenium` and `10x_visium_hd`, `--in-json` can replace these flags (see [Parameters](#sge-conversion-parameters)).

    <sub><sup>1</sup>: For SGE from platforms not explicitly supported, or from custom/preprocessed sources. The default columns are `X`, `Y`, `gene` in a tab-delimited file; override them with the `--csv-*` options.

!!! info "Pre-installed tools"
    * `gzip` (or `pigz`).
    * `spatula`: required for every MEX conversion (`10x_visium_hd`, `illumina`, and `seqscope` with `--in-mex`) and for `--sge-visual`. Defaults to the repository submodule (`<cartloader_dir>/submodules/spatula/bin/spatula`).
    * Parquet input (`--in-parquet`, `--pos-parquet`): converted with `cartloader parquet_to_csv_rapid` (Python `polars` + `pigz`), or with `parquet-tools` when `--use-parquet-tools` is set.
    * `--filter-by-density`: the `ficture` Python package (installed with CartLoader), which provides `ficture filter_by_density`.
    * `--north-up`: `gdal_translate` and `gdalwarp`.

---

## Actions

!!! warning "Action Specifications"
    SGE conversion always runs. `--filter-by-density` and `--sge-visual` are optional add-ons.

### SGE Conversion (always runs)
Converts the input into a FICTURE-compatible transcript TSV, rescaling coordinates to µm (`--units-per-um`). Feature filters (lists or regexes) and, for 10x Xenium, a Phred quality filter are applied during conversion.

### Density-based Filtering (`--filter-by-density`)
Identifies tissue regions from transcript density and spatial structure (`ficture filter_by_density`) and writes a filtered SGE alongside the unfiltered one. The filtered SGE then becomes the SGE recorded in `sge_assets.json`.

### SGE Visualization (`--sge-visual`)
Draws a monochrome image of the transcript positions with `spatula draw-xy`. Add `--north-up` to also write a georeferenced, north-up GeoTIFF (Y increases northward, X eastward). With `--filter-by-density`, both the unfiltered and the filtered SGE are drawn.

---
## Example Usage

### 1) Input SGE in MEX Format

#### 1.1) Seq-Scope
```bash
cartloader sge_convert \
    --platform seqscope \
    --in-mex /path/to/input/mex/dir \
    --units-per-um 1000 \
    --icols-mtx 1 \
    --colname-count count \
    --out-dir /path/to/output/dir \
    --filter-by-density \
    --sge-visual \
    --north-up
```

For the raw per-molecule Seq-Scope TSV, pass `--in-csv /path/to/raw.tsv.gz` instead of `--in-mex` (no `--units-per-um` needed: the TSV route defaults to `1000`).

#### 1.2) 10x Visium HD
```bash
cartloader sge_convert \
    --platform 10x_visium_hd \
    --in-mex /path/to/binned_outputs/square_002um/filtered_feature_bc_matrix \
    --pos-parquet /path/to/binned_outputs/square_002um/spatial/tissue_positions.parquet \
    --scale-json /path/to/binned_outputs/square_002um/spatial/scalefactors_json.json \
    --exclude-feature-regex '^(BLANK|NegCon|NegPrb)' \
    --out-dir /path/to/output/dir \
    --sge-visual
```

#### 1.3) Illumina StrataMap
```bash
cartloader sge_convert \
    --platform illumina \
    --in-mex /path/to/input/mex/dir \
    --out-dir /path/to/output/dir
```

The barcode delimiter (`:`), the coordinate fields (x = field 3, y = field 2) and `--units-per-um` (nanometer vs. micron barcodes) are set automatically.

### 2) Input SGE in TSV/CSV Format

This applies to every CSV/TSV platform. Below converts a Stereo-seq GEM file (500 nm bins, i.e. 2 units per µm; `--csv-comment` skips the `#` metadata lines at the top of the file):

```bash
cartloader sge_convert \
    --platform bgi_stereoseq \
    --in-csv /path/to/input.gem.gz \
    --csv-comment \
    --units-per-um 2 \
    --out-dir /path/to/output/dir \
    --exclude-feature-regex '^(mt-|Gm[0-9])' \
    --filter-by-density \
    --sge-visual \
    --north-up
```

!!! warning "Verify Your Input Structure"
    Always check the column names and delimiter of your input. If they differ from the platform defaults below, override them with the `--csv-*` options (and `--min-phred-score` for quality filtering).

!!! info "Platform-specific Default Settings (CSV/TSV input)"

    | Platform                     | `--csv-comment`<sup>1</sup> | `--csv-delim` | `--csv-colname-x` | `--csv-colname-y` | `--csv-colname-count` | `--csv-colname-feature-name` |
    |------------------------------|-----------------------------|---------------|-------------------|-------------------|-----------------------|------------------------------|
    | 10x Xenium<sup>2</sup>       | `False`                     | `,`           | `x_location`      | `y_location`      | -                     | `feature_name`               |
    | Stereo-seq                   | `False`                     | `\t`          | `x`               | `y`               | `ExonCount`           | `geneName`                   |
    | CosMx SMI<sup>3</sup>        | `False`                     | `,`           | `x_global_px`     | `y_global_px`     | -                     | `target`                     |
    | MERSCOPE                     | `False`                     | `,`           | `global_x`        | `global_y`        | -                     | `gene`                       |
    | Pixel-seq                    | `False`                     | `\t`          | `xcoord`          | `ycoord`          | -                     | `geneName`                   |
    | Nova-ST                      | `True`                      | `\t`          | `x`               | `y`               | `MIDCount`            | `geneID`                     |
    | Seq-Scope (`--in-csv`)<sup>4</sup> | `False`               | `\t`          | `X`               | `Y`               | `gn`                  | `gene`                       |
    | Generic                      | `False`                     | `\t`          | `X`               | `Y`               | -                     | `gene`                       |

    `-`: no count column; each row counts as 1. With a count column, rows whose count is 0 are dropped.

    <sub><sup>1</sup> `--csv-comment`: if `True`, lines starting with `#` are skipped. Without it, a `#`-prefixed header line (e.g. the Seq-Scope TSV's `#lane ...`) is read as the header with the `#` stripped.</sub>

    <sub><sup>2</sup> 10x Xenium also defaults to `--csv-colname-phredscore qv`, `--min-phred-score 20`, and `--exclude-feature-regex '^(BLANK|Blank-|Deprecated|Intergenic|NegCon|Unassigned)'`.</sub>

    <sub><sup>3</sup> CosMx global coordinates are in pixels, so set `--units-per-um` to pixels per µm. `run_together` ingests CosMx with [`reformat_cosmx`](./reformat_cosmx.md) instead.</sub>

    <sub><sup>4</sup> The Seq-Scope TSV route also defaults to `--units-per-um 1000` and `--csv-colnames-positive lane tile X Y`.</sub>

---
## Parameters

Below are the core parameters. See more details in the collapsible sections below.

### SGE Conversion Parameters

* `--platform` (str, required): Source platform. Options: `10x_visium_hd`, `seqscope`, `illumina`, `10x_xenium`, `bgi_stereoseq`, `cosmx_smi`, `vizgen_merscope`, `pixel_seq`, `nova_st`, `generic`.
* `--out-dir` (str, required): Output directory.
* `--in-mex` (str): Input MEX directory (platforms: `10x_visium_hd`, `seqscope`, `illumina`; default: current working directory).
* `--in-csv` (str): Input CSV/TSV (CSV/TSV platforms; for `seqscope` it selects the raw per-molecule TSV route instead of `--in-mex`).
* `--in-parquet` (str): Input transcript parquet (platform: `10x_xenium`).
* `--pos-parquet` (str): Position parquet with spatial coordinates (platform: `10x_visium_hd`, required; typical: `tissue_positions.parquet`).
* `--scale-json` (str): Scale JSON; sets `--units-per-um` to `1 / microns_per_pixel` (platform: `10x_visium_hd`; typical: `scalefactors_json.json`).
* `--in-json` (str): Input manifest JSON whose `SGE` block names the inputs (as written by `load_xenium_ranger` / `load_space_ranger`). Replaces `--in-csv`/`--in-parquet` for `10x_xenium`, and `--in-mex`/`--pos-parquet`/`--scale-json` for `10x_visium_hd` (where it also sets `--scale-xy 2`).
* `--units-per-um` (float): Coordinate units per µm in the input (default: `1`). Platform defaults: `1000` for the Seq-Scope TSV route; for `illumina` it is **detected from the barcode file** (`1000` for nanometer barcodes such as `SBC:433503:2393851`, `1` for micron barcodes such as `SBC:686.951:4668.15`), and an explicit value overrides the detection (a contradicting value is reported as a warning). For `10x_visium_hd`, `--scale-json` takes precedence over this flag.

### Feature Filtering Parameters

* `--include-feature-regex` (regex): Keep only features whose name matches.
* `--exclude-feature-regex` (regex): Drop features whose name matches (default for `10x_xenium`: `'^(BLANK|Blank-|Deprecated|Intergenic|NegCon|Unassigned)'`; none for other platforms). Pass `''` to disable the default.
* `--include-feature-list` (str): File listing the feature names to keep (one per line).
* `--exclude-feature-list` (str): File listing the feature names to drop (one per line).

Regexes use **substring (search)** semantics, so anchor them with `^`/`$` as needed. On the MEX platforms (`10x_visium_hd`, `illumina`, and `seqscope` with `--in-mex`), filtering runs inside `spatula convert-sge`, which accepts at most one include-type and one exclude-type filter: combining a list and a regex of the same polarity is an error there. Merge them into one list with [`cartloader feature_select`](./feature_customization.md#resolving-filters-into-one-list-feature_select). On the CSV/TSV platforms all four filters combine.

??? note "Auxiliary SGE Conversion Parameters"

    Recommend to use the default values; override only if needed. See more details by running:
    ```bash
    cartloader sge_convert --help
    ```

    **Auxiliary Input MEX Parameters** (`10x_visium_hd`, `seqscope`, `illumina`):

    * `--mex-bcd`, `--mex-ftr`, `--mex-mtx` (str): File names of barcodes, features and matrix inside `--in-mex` (defaults: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`).
    * `--icols-mtx` (int): 1-based index of the count column to use in the matrix file (default: 1).
    * `--icol-ftr-id`, `--icol-ftr-name` (int): 1-based columns of feature ID and name in the feature file (defaults: 1, 2).
    * `--icol-bcd-barcode` (int): 1-based column of the barcode in the barcode file (platform: `seqscope`; default: 1).
    * `--icol-bcd-x`, `--icol-bcd-y` (int): 1-based column (Seq-Scope) or barcode field (Illumina) holding X and Y (defaults: 6 and 7 for `seqscope`; 3 and 2 for `illumina`).
    * `--bcd-delim` (str): Delimiter splitting an Illumina barcode into fields (platform: `illumina`; default: `:`).
    * `--pos-colname-barcode`, `--pos-colname-x`, `--pos-colname-y` (str): Columns of `--pos-parquet` (platform: `10x_visium_hd`; defaults: `barcode`, `pxl_col_in_fullres`, `pxl_row_in_fullres`).
    * `--pos-delim` (str): Delimiter of the converted position file (platform: `10x_visium_hd`; default: `,`).
    * `--allow-duplicate-gene-names` (flag): Allow duplicate gene names (platform: `10x_visium_hd`).
    * `--keep-mismatches` (flag): Keep barcodes that did not match a spatial position (X = Y = 0); by default they are dropped by `cartloader sge_drop_mismatches` (platform: `seqscope` with `--in-mex`; for debugging).

    **Auxiliary Input CSV/TSV Parameters** (defaults per platform: see the table above):

    * `--csv-comment` (flag): Skip lines starting with `#`.
    * `--csv-delim` (str): Delimiter of the input file.
    * `--csv-colname-x`, `--csv-colname-y` (str): Columns of X and Y coordinates.
    * `--csv-colname-count` (str): Column of the expression count. Without it, each row counts as 1.
    * `--csv-colname-feature-name` (str): Column of the gene name.
    * `--csv-colnames-others` (list): Additional columns to keep in the output (e.g., `cell_id overlaps_nucleus`).
    * `--csv-colnames-positive` (list): Columns that must be strictly positive; rows with a zero, negative or non-numeric value are dropped as malformed (default: `lane tile X Y` for the Seq-Scope TSV route).
    * `--csv-colname-phredscore` (str): Column of the Phred-scaled quality score (default: `qv` for `10x_xenium`).
    * `--min-phred-score` (float): Minimum quality score to keep a transcript (default: `20` for `10x_xenium`; no filtering for others).

    **Auxiliary Output Parameters**:

    * `--out-transcript` (str): File name of the output transcript TSV (default: `transcripts.unsorted.tsv.gz`).
    * `--out-minmax` (str): File name of the coordinate min/max TSV (default: `coordinate_minmax.tsv`).
    * `--out-feature` (str): File name of the per-gene count TSV (default: `feature.clean.tsv.gz`).
    * `--out-json` (str): Output JSON manifest of SGE paths (default: `<out-dir>/sge_assets.json`).
    * `--precision-um` (int): Decimal places for output coordinates; `0` rounds to integers (default: 2).
    * `--colname-x`, `--colname-y` (str): Output X/Y column names (defaults: `X`, `Y`).
    * `--colname-count` (str): Output count column name (default: `count`).
    * `--colname-feature-name` (str): Output gene column name (default: `gene`).
    * `--jitter-xy` (float): Maximum random jitter (µm) added to each coordinate (default: 0, off).

    **Environment Parameters**:

    * `--gzip` (str): Path to `gzip` (default: `gzip`; e.g., `"pigz -p 4"` for speed).
    * `--spatula` (str): Path to `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`).
    * `--pigz` (str), `--pigz-threads` (int): `pigz` binary and threads used for parquet conversion (defaults: `pigz`, 4).
    * `--use-parquet-tools` (flag): Convert parquet with `parquet-tools` instead of `polars`/`pigz` (may be slower for large files).
    * `--parquet-tools` (str): Path to `parquet-tools` (used only with `--use-parquet-tools`; default: `parquet-tools`).

    **Run Parameters**:

    * `--dry-run` (flag): Write the Makefile and print the commands without executing.
    * `--restart` (flag): Ignore existing outputs and rerun all steps.
    * `--makefn` (str): Makefile name under `--out-dir` (default: `sge_convert.mk`).
    * `--n-jobs`, `-j` (int): Number of parallel jobs (default: 1).

### Density‑based Filtering Parameters

* `--filter-by-density` (flag): Enable density-based filtering.
* `--out-filtered-prefix` (str): Prefix of the filtered outputs under `--out-dir` (default: `filtered`).

??? note "Auxiliary Density‑based Filtering Parameters"
    We recommend using default values; override only if needed.

    * `--radius` (int): Radius for the polygon area calculation (default: 15).
    * `--quartile` (int): Quartile for the polygon area calculation (default: 2).
    * `--hex-n-move` (int): Sliding step (default: 1).
    * `--polygon-min-size` (int): Minimum polygon size (default: 500).

### SGE Visualization Parameters
* `--sge-visual` (flag): Draw the SGE image.
* `--north-up` (flag): Also write a north-up GeoTIFF (with `--sge-visual`).

??? note "Auxiliary SGE Visualization Parameters"
    We recommend using default values; override only if needed.

    * `--out-xy` (str): File name of the SGE image (default: `xy.png`; the filtered image is prefixed by `--out-filtered-prefix`).
    * `--out-northup-tif` (str): File name of the north-up image (default: `xy_northup.tif`).
    * `--scale-xy` (float): µm per pixel of the image (default: 1.0).
    * `--srs` (str): Spatial reference system (with `--north-up`; default: `EPSG:3857`).
    * `--resample` (str): Resampling method (with `--north-up`; e.g., `near`, `bilinear`, `cubic`; default: `cubic`).
    * `--gdal_translate`, `--gdalwarp` (str): Paths to the GDAL binaries (with `--north-up`; defaults: `gdal_translate`, `gdalwarp`).

---

## Output

All outputs are written under `--out-dir`.

### Unified SGE matrix

SGE conversion writes:

{%
  include-markdown "../../includes/includemd_vigenettes_sgeformat.md"
%}

With `--filter-by-density`, the filtered SGE uses fixed names (independent of `--out-transcript`/`--out-feature`/`--out-minmax`):

* `<prefix>.transcripts.unsorted.tsv.gz`, `<prefix>.feature.lenient.tsv.gz`, `<prefix>.coordinate_minmax.tsv` (`<prefix>` = `--out-filtered-prefix`)
* tissue boundaries from `ficture filter_by_density` (`<prefix>.boundary.lenient.geojson`, `<prefix>.boundary.strict.geojson`, and PNG previews)

### SGE manifest (`sge_assets.json`)

`sge_convert` records the SGE to use downstream in `sge_assets.json` (or `--out-json`), written when the Makefile is generated:

```json
{
  "transcript": "/path/to/out/transcripts.unsorted.tsv.gz",
  "feature": "/path/to/out/feature.clean.tsv.gz",
  "minmax": "/path/to/out/coordinate_minmax.tsv",
  "density_filtering": false,
  "visual": { "png": "/path/to/out/xy.png", "northup": false }
}
```

* `transcript`, `feature`, `minmax`: the SGE files. With `--filter-by-density` they point to the filtered files and `density_filtering` is `true`.
* `visual` (only with `--sge-visual`): the image (`png` is the north-up TIFF when `--north-up` is set, and the filtered image when filtering is enabled).

[`run_ficture2 --in-json`](./run_ficture2.md) reads this file directly.

### SGE Images

* `--sge-visual`: monochrome PNG of the SGE (`xy.png`, plus `<prefix>.xy.png` when filtered).
* `--north-up`: georeferenced north-up GeoTIFF (`xy_northup.tif`, plus `<prefix>.xy_northup.tif` when filtered).

### Other files

* `sge_convert.done` (and `sge_density_filtering.done`): done flags used to skip completed steps on rerun.
* `sge_convert.mk`: the generated Makefile.
* `transcripts.parquet.csv.gz`: CSV converted from `--in-parquet` (10x Xenium only).
