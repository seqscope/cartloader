# Spatial Gene Expression Format Conversion

## Overview

SGE datasets vary widely in format and resolution across platforms. Since `FICTURE` requires SGE in a specific format, `CartLoader` toolkit provides the `sge_convert` module to standardize raw, platform-specific SGE into a FICTURE-compatible format.

## Requirements

!!! info "Input Data Requirements"
    Ensure the input data (raw, platform‑specific SGE) is transcript‑indexed and contains at least the following fields:

    * Spatial coordinates (X, Y)
    * Feature metadata (such as gene symbols)
    * Expression counts

!!! info "Platform Compatibility"

    The current `sge_convert` supports standardizing SGE from the following platforms:

    | Source                                                                            | `--platform` Option | Required Input Files                                |
    |-----------------------------------------------------------------------------------|---------------------|-----------------------------------------------------|
    | [10x Visium HD](https://www.10xgenomics.com/platforms/visium)                     | `10x_visium_hd`     | `--in-mex`, `--pos-parquet`, `--scale-json`         |
    | [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0)                   | `seqscope`          | `--in-mex`                                          |
    | [10x Xenium](https://www.10xgenomics.com/platforms/xenium)                        | `10x_xenium`        | `--in-csv` (common); `--in-parquet` also supported. |
    | [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq) | `bgi_stereoseq`     | `--in-csv`                                          |
    | [CosMx SMI](https://nanostring.com/products/cosmx-spatial-molecular-imager)       | `cosmx_smi`         | `--in-csv`                                          |
    | [Vizgen MERSCOPE](https://vizgen.com/merscope-ultra/)                             | `vizgen_merscope`   | `--in-csv`                                          |
    | [Pixel-seq](https://www.cell.com/cell/fulltext/S0092-8674(22)01367-8)             | `pixel_seq`         | `--in-csv`                                          |
    | Nova-ST                                                                           | `nova_st`           | `--in-csv`                                          |
    | Generic CSV/TSV input <sup>1</sup>                                                | `generic`           | `--in-csv`                                          |

    <sub><sup>1</sup>: For SGE from platforms not yet explicitly supported by `CartLoader`, or from custom/preprocessed sources, `sge_convert` provides a `generic` option that accepts CSV/TSV files with basic required fields (e.g., gene, spatial coordinates, expression count) for standardization and processing.

!!! info "Pre-installed tools"
    * `gzip` (or pigz)
    * `spatula` (required if `--sge-visual` is set)
    * `gdal_translate`, `gdalwarp` (required if `--sge-visual` is set with `--north-up`)


## Example Usage

### 1) Input SGE in MEX Format

#### 1.1) `Seq-Scope`
```bash
cartloader sge_convert \
    --platform seqscope \
    --in-mex /path/to/input/dir/of/mex \   
    --units-per-um 1000 \
    --icols-mtx 1 \
    --out-dir /path/to/output/dir \
    --colnames-count count \
    --filter-by-density \
    --out-filtered-prefix filtered \
    --genomic-feature count \
    --sge-visual \
    --north-up \
    --spatula /path/to/spatula/binary 
```

#### 1.2) `10X Visium HD`
```bash
cartloader sge_convert \
    --platform 10x_visium_hd \
    --in-mex /path/to/input/dir/of/mex \   
    --in-parquet /path/to/input/parquet/file \
    --scale-json /path/to/input/json/file \
    --exclude-feature-regex '^(BLANK.*$|NegCon.*$|NegPrb.*$)' \
    --out-dir /path/to/output/dir \
    --spatula /path/to/spatula/binary \
    --sge-visual 
```

### 2) Input SGE in TSV/CSV Format

This applies to input SGE in TSV/CSV format from platforms including 10x Xenium, Stereo‑seq, CosMx SMI, MERSCOPE, Pixel‑seq, and Nova‑ST. To simplify preprocessing, `sge_convert` automatically applies platform‑specific defaults for common CSV/TSV parameters.

Below is an example converting SGE from Stere-seq.

```bash
cartloader sge_convert \
    --platform bgi_stereoseq \
    --in-csv /path/to/input/csv/file \
    --units-per-um 28.75 \
    --out-dir /path/to/output/dir \
    --colnames-count count \
    --exclude-feature-regex '^(BLANK.*$|NegCon.*$|NegPrb.*$)' \
    --filter-by-density \
    --out-filtered-prefix filtered \
    --genomic-feature count \
    --sge-visual \
    --north-up \
    --spatula /path/to/spatula/binary 
```

!!! warning "Verify Your Input Structure"
    Always double-check the column names and metadata format of your input files. If they differ from the expected defaults, override them using `--csv-*` and `--min-phred-score` options.

??? "Click to view platform-specific default settings"

    To streamline the process, `sge_convert` automatically applies platform‑dependent defaults for CSV/TSV parsing based on known file formats and column conventions.
    Below summarizes the default values for key parameters per supported platform:

    | Platform               | `--csv-comment`<sup>1</sup> | `--csv-delim` | `--csv-colname-x` | `--csv-colname-y` | `--csv-colnames-count` | `--csv-colname-feature-name` |
    |------------------------|-----------------------------|---------------|-------------------|-------------------|------------------------|------------------------------|
    | 10X Xenium<sup>2</sup> | `False`                     | `,`           | `x_location`      | `y_location`      | -                      | `feature_name`               |
    | StereoSeq              | `False`                     | `\t`          | `x`               | `y`               | `MIDCounts`            | `geneID`                     |
    | CosMx SMI              | `False`                     | `,`           | `x_local_px`      | `y_local_px`      | -                      | `target`                     |
    | MERSCOPE               | `False`                     | `,`           | `global_x`        | `global_y`        | -                      | `gene`                       |
    | Pixel-seq              | `False`                     | `\t`          | `xcoord`          | `ycoord`          | -                      | `geneName`                   |
    
    <sub><sup>1</sup> `--csv-comment`: If `True`, the lines starts with `#` will be treated as comments and will be skipped.

    <sub><sup>2</sup> 10X Xenium: Besides the above default settings, for 10X Xenium data, `sge_convert` also applies `--csv-colname-phredscore qv ` and `--min-phred-score 20`.


---
## Actions

### SGE Conversion
Converting SGE into a FICTURE-compatible TSV format. During conversion, SGE coordinates are rescaled to micrometer units based on the pixel resolution specified in the input. It's also available to apply feature (typically, genes) filtering.

### Density-based Filtering (`--filter-by-density`)
If `--filter-by-density` is set, automatically identify and retain high-quality tissue regions based on transcript density and spatial structure. This step takes the format-standardized SGE as input and generate a density-based filtered SGE.

### SGE Visualization (`--sge-visual`)
If `--sge-visual` is set, draws an image of 2D points provided as an input. In this step, it is optional to enable the `--north-up` option to ensuring correct spatial orientation (i.e., Y-axis increases upward/north and X-axis increases to the right/east).

---
## Parameters

Below are the core parameters. See more details in the collapsible sections below.

### SGE Conversion

* `--platform` (str): Source platform to infer input format and defaults. Options: `10x_visium_hd`, `seqscope`, `10x_xenium`, `bgi_stereoseq`, `cosmx_smi`, `vizgen_merscope`, `pixel_seq`, `nova_st`, `generic`.
* `--in-json` (str): Input manifest JSON. If set, can skip `--in-mex/--in-parquet/--in-csv/--pos-parquet/--scale-json` (platforms: 10x_xenium, 10x_visium_hd).
* `--in-mex` (str): Path to input MEX directory (platforms: 10x Visium HD, SeqScope).
* `--in-csv` (str): Path to input CSV/TSV (platforms: 10x Xenium, BGI Stereo‑seq, CosMx SMI, Vizgen MERSCOPE, Pixel‑seq, Nova‑ST, generic).
* `--in-parquet` (str): Path to input transcript parquet (platform: 10x Xenium).
* `--pos-parquet` (str): Path to position parquet with spatial coordinates (platform: 10x Visium HD; typical: `tissue_positions.parquet`).
* `--scale-json` (str): Path to scale JSON; if set, derives `--units-per-um` from `microns_per_pixel` (platform: 10x Visium HD; typical: `scalefactors_json.json`).
* `--units-per-um` (float): Coordinate units per µm (default: 1.00). Prefer `--scale-json` for 10x Visium HD.
* `--out-dir` (str): Output directory.
* `--include-feature-regex` (regex): Regex of feature/gene names to include.
* `--exclude-feature-regex` (regex): Regex of feature/gene names to exclude.
<!-- * `--include-feature-list` and `--exclude-feature-list`: A file containing a list of input genes to be included or excluded. 
* `--include-feature-type-regex`: A regex pattern of feature/gene type to be included.
* `--csv-colname-feature-type`: If `--include-feature-type-regex`, column name in the input TSV/CSV that contains gene type information.
* `--feature-type-ref`, `--feature-type-ref-delim`, `--feature-type-ref-colidx-name`, and `--feature-type-ref-colidx-type`: If `--include-feature-type-regex`, define the file path, delimiter, and column indices for gene names and types of an additional reference file to provide gene type information. -->
??? note "Auxiliary SGE Conversion Paramaters"

    Recommend to use the default values; override only if needed. See more details by running:
    ```bash
    CartLoadersge_convert --help
    ```

    **Auxiliary Input MEX Parameters**:

    * `--icols-mtx` (int or comma-spearated list): Comma-separated, 1-based indices of the target genomic features among the count columns in the input matrix file (default: 1)
    * `--colnames-count` (string or comma-spearated list): Comma-separated output column names for the specified genomic features (default: count). The number of names specified by `--colnames-count` must match the number of indices provided in `--icols-mtx`.

    **Auxiliary Input CSV/TSV Parameters**:

    * `--csv-comment` (flag): If enabled, lines starting with `#` are skipped (default: `False` for 10x Xenium, Stereo‑seq, CosMx SMI, MERSCOPE, and Pixel‑seq).
    * `--csv-delim` (str): Delimiter for the input file (default: `","` for 10x Xenium, CosMx SMI, and MERSCOPE; `"\t"` for Stereo‑seq, Pixel‑seq).
    * `--csv-colname-x` (str): Column name for X coordinates (default: `x_location` for 10x Xenium; `x` for Stereo‑seq; `x_local_px` for CosMx SMI; `global_x` for MERSCOPE; `xcoord` for Pixel‑seq).
    * `--csv-colname-y` (str): Column name for Y coordinates (default: `y_location` for 10x Xenium; `y` for Stereo‑seq; `y_local_px` for CosMx SMI; `global_y` for MERSCOPE; `ycoord` for Pixel‑seq).
    * `--csv-colnames-count` (str): Comma‑separated column names for expression count. If not provided, defaults to a count of 1 per transcript (default: `MIDCounts` for Stereo‑seq).
    * `--csv-colname-feature-name` (str): Column name for gene name (default: `feature_name` for 10x Xenium; `geneID` for Stereo‑seq; `target` for CosMx SMI; `gene` for MERSCOPE; `geneName` for Pixel‑seq).
    * `--csv-colnames-others` (str): Column names to keep.
    * `--csv-colname-phredscore` (str): Column name for Phred‑scaled quality value estimating the probability of incorrect calls (default: `qv` for 10x Xenium).
    * `--min-phred-score` (int): Phred‑scaled quality score cutoff (default: `20` for 10x Xenium).

    **Auxiliary Output Parameters**: 

    * `--out-transcript` (str): File name for output compressed transcript-indexed SGE file in TSV format (default: `transcripts.unsorted.tsv.gz`).
    * `--out-minmax` (str): File name for coordinate min–max values in TSV format (default: `coordinate_minmax.tsv`).
    * `--out-feature` (str): File name for compressed UMI count per gene in TSV format (default: `feature.clean.tsv.gz`).
    * `--precision-um` (int): Decimal precision for transcript coordinates; set to `0` to round to integers (default: 2).
    * `--colname-x` (str): Column name for the X-coordinate in the output SGE (default: X).
    * `--colname-y` (str): Column name for the Y-coordinate in the output SGE (default: Y).
    * `--colname-count` (str): Comma‑separated column names for count in the output SGE (default: count).
    * `--colname-feature-name` (str): Column name for the gene name in the output SGE (default: gene).
    * `--out-json` (str): Output JSON manifest of SGE paths (default: `<out-dir>/sge_assets.json`).

    **Environment Parameters**  
    If the binaries are already available in your system's `PATH`, you may omit these options.

    * `--gzip` (str): Path to `gzip` binary; consider `pigz -p 4` for faster processing (default: `gzip`)
    * `--spatula` (str): Path to `spatula` binary (default: `spatula`)
    * `--parquet-tools` (str): Path to `parquet-tools` binary (used with --in-parquet or --pos-parquet; default: `parquet-tools`)

    **Run Parameters**:

    * `--dry-run` (flag): Generate the Makefile but do not execute it.
    * `--restart` (flag): Ignore existing outputs and re-run all steps.
    * `--makefn` (str): Output Makefile name (default: `sge_convert.mk`)
    * `--n-jobs` (int): Parallel jobs when executing the Makefile (default: 1).

### Density‑based Filtering

* `--filter-by-density` (flag): Enable SGE filtering by density.
* `--out-filtered-prefix` (str): Prefix for output filtered SGE files (default: filtered).

??? note "Auxiliary Density‑based Filtering Parameters"
    We recommend using default values; override only if needed.

    * `--radius` (int): Radius for the polygon area calculation (default: 15).
    * `--quartile` (int): Quartile for the polygon area calculation (default: 2).
    * `--hex-n-move` (int): Sliding step (default: 1).
    * `--polygon-min-size` (int): Minimum polygon size (default: 500).

### SGE Visualization
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

---

## Output

`CartLoader` generates the following harmonized outputs:

### Unified SGE matrix

Both SGE conversion and density-based filtering generate a unified SGE matrix, consisting of:

{%
  include-markdown "../../includes/includemd_vigenettes_sgeformat.md"
%}

### SGE Images

* When `--sge-visual` is enabled, a monochrome PNG image is generated to visualize the SGE data.
* When `--north-up` is enabled, a georeferenced TIFF image is produced with a north-up orientation.
