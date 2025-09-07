# Spatial Gene Expression Format Conversion

## Overview

SGE datasets vary widely in format and resolution across platforms. Since `FICTURE` requires SGE in a specific format, `cartloader` toolkit provides this `sge_convert` module to standardize the raw, platform-specific SGE into a FICTURE-compatible format.

## Requirements

!!! info "Input Data Requirement:"
    Please make sure the input data (raw, platform-specific SGE) contains at least the following information. should be a transcript-indexed SGE containing at least:

    * Spatial coordinates (X coordinates, Y coordinates)
    * Feature metadata (such as gene symbols)
    * Expression Counts

!!! info "Platform Compatibility"

    The current `sge_convert` supports standarizing SGE from the following platforms:

    | Source                                                                            | `--platform` Option | Required Input Files                       |
    |-----------------------------------------------------------------------------------|---------------------|--------------------------------------------|
    | [10x Visium HD](https://www.10xgenomics.com/platforms/visium)                     | `10x_visium_hd`     | `--in-mex`, `--in-parquet`, `--scale-json` |
    | [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0)                   | `seqscope`          | `--in-mex`                                 |
    | [10x Xenium](https://www.10xgenomics.com/platforms/xenium)                        | `10x_xenium`        | `--in-csv`                                 |
    | [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq) | `bgi_stereoseq`     | `--in-csv`                                 |
    | [CosMx SMI](https://nanostring.com/products/cosmx-spatial-molecular-imager)       | `cosmx_smi`         | `--in-csv`                                 |
    | [Vizgen MERSCOPE](https://vizgen.com/merscope-ultra/)                             | `vizgen_merscope`   | `--in-csv`                                 |
    | [Pixel-seq](https://www.cell.com/cell/fulltext/S0092-8674(22)01367-8)             | `pixel_seq`         | `--in-csv`                                 |
    | Generic CSV/TSV input <sup>1</sup>                                                | `generic`           | `--in-csv`                                 |

    <sub><sup>1</sup>: For SGE from platforms not yet explicitly supported by `cartloader`, or from custom/preprocessed sources, `sge_convert` provides a `generic` option that accepts CSV/TSV files with basic required fields (e.g., gene, spatial coordinates, expression count) for standardization and processing.

## Example Usages

### Input SGE in MEX Format

#### `Seq-Scope`
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
    --spatula /path/to/spatul/binary 
```

#### `10X Visium HD`
```bash
cartloader sge_convert \
    --platform 10x_visium_hd \
    --in-mex /path/to/input/dir/of/mex \   
    --in-parquet /path/to/input/parquet/file \
    --scale-json /path/to/input/json/file \
    --exclude-feature-regex '^(BLANK.*$|NegCon.*$|NegPrb.*$)' \
    --out-dir /path/to/output/dir \
    --spatula /path/to/spatul/binary \
    --sge-visual 
```

### Input SGE in TSV/CSV Format

This applies to input SGE in TSV/CSV format from platforms including 10X Xenium, StereoSeq, Cosmx SMI, MERSCOPE, Pixel-seq. To simplify preprocessing, `sge_convert` automatically applies platform-specific defaults for common CSV/TSV parameters.

Below is an example converting SGE from StereoSeq.

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
    --spatula /path/to/spatul/binary 
```

!!! warning "Verify your input structure"
    Always double-check the column names and metadata format of your input files. If they differ from the expected defaults, override them using `--csv-*` and `--min-phred-score` options.

??? "Click to view platform-specific default settings"

    To streamline the process, `sge_convert` automatically applies platform-dependent defaults for CSV/TSV parsing based on known file formats and column conventions.
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

### (Optional) Density-based Filtering
Automatically identify and retain high-quality tissue regions based on transcript density and spatial structure. This step takes the format-standardized SGE as input and generate a density-based filtered SGE.

### (Optional) SGE Visualization
Draws an image of 2D points provided as an input. In this step, it is optional to enable the `--north-up` option to ensuring correct spatial orientation (i.e., Y-axis increases upward/north and X-axis increases to the right/east).

---
## Parameters

The following outlines the **minimum required parameters**.

For most auxiliary parameters, the **default values are recommended** and could be modified when they do not suit your use case. See more details in the collapsible sections below or by running:
```bash
cartloader sge_convert --help
```

### SGE Conversion

* `--platform` (str): Source platform to infer the input file format and default setting (options: "`10x_visium_hd`", "`seqscope`", "`10x_xenium`", "`bgi_stereoseq`", "`cosmx_smi`", "`vizgen_merscope`", "`pixel_seq`", "`generic`").
* `--in-mex` (str): Path to the input SGE directory in MEX format. Required for MEX-formatted data (e.g., 10X Visium HD, SeqScope).
* `--in-csv` (str): Path to the input SGE file in CSV or TSV format. Required for CSV/TSV-formatted data (e.g., 10X Xenium, StereoSeq, CosMX SMI, MERSCOPE, Pixel-seq).
* `--in-parquet` (str) : Path to the input parquet file with spatial coordinates, if available. Typically named `tissue_positions.parquet` (10X Visium HD).
* `--scale-json` (str): Path to the scale JSON file to compute `--units-per-um`, if available. Typically named `scalefactors_json.json` (10X Visium HD).
* `--units-per-um` (float): Coordinate units per micrometer (default: 1.00). Skip if `--scale-json` is provided.
* `--out-dir` (str) : Path to the output directory.
* `--include-feature-regex` (regex): (Optional) A regex pattern of feature/gene names to be included.
* `--exclude-feature-regex` (regex): (Optional) A regex pattern of feature/gene names to be excluded.
<!-- * `--include-feature-list` and `--exclude-feature-list`: A file containing a list of input genes to be included or excluded. 
* `--include-feature-type-regex`: A regex pattern of feature/gene type to be included.
* `--csv-colname-feature-type`: If `--include-feature-type-regex`, column name in the input TSV/CSV that contains gene type information.
* `--feature-type-ref`, `--feature-type-ref-delim`, `--feature-type-ref-colidx-name`, and `--feature-type-ref-colidx-type`: If `--include-feature-type-regex`, define the file path, delimiter, and column indices for gene names and types of an additional reference file to provide gene type information. -->
??? note "Auxiliary SGE Conversion Paramaters"

    **Auxiliary Input MEX Parameters**:

    * `--icols-mtx` (int or comma-spearated list): Comma-separated, 1-based indices of the target genomic features among the count columns in the input matrix file. (Default: 1)
    * `--colnames-count` (string or comma-spearated list): Comma-separated output column names for the specified genomic features. (Default: count). The number of names specified by `--colnames-count` must match the number of indices provided in `--icols-mtx`.

    **Auxiliary Input CSV/TSV Parameters**:

    * `--csv-comment` (flag): If enabled, lines starts with `#` will be skipped (default: `False` for 10X Xenium, StereoSeq, CosMx SMI, MERSCOPE, and Pixel-seq).
    * `--csv-delim` (str): Delimiter for the input file (default: `","` for 10X Xenium, CosMx SMI, and MERSCOPE; `"\t"` for StereoSeq, Pixel-seq).
    * `--csv-colname-x` (str): Column name for X coordinates (default: `x_location` for 10X Xenium; `x` for StereoSeq; `x_local_px` for CosMx SMI; `global_x` for MERSCOPE; `xcoord` for Pixel-seq).
    * `--csv-colname-y` (str): Column name for Y coordinates (default: `y_location` for 10X Xenium; `y` for StereoSeq; `y_local_px` for CosMx SMI; `global_y` for MERSCOPE; `ycoord` for Pixel-seq).
    * `--csv-colnames-count` (str): Comma-separated column names for expression count. If not provided, a count of 1 per transcript (default: `MIDCounts` for StereoSeq).
    * `--csv-colname-feature-name` (str): Column name for gene name (default: `feature_name` for 10X Xenium; `geneID` for StereoSeq; `target` for CosMx SMI; `gene` for MERSCOPE; `geneName` for Pixel-seq).
    * `--csv-colnames-others` (str): Columns names to keep.
    * `--csv-colname-phredscore` (str): Column name for Phred-scaled quality value estimating the probability of incorrect calls (default: `qv` for 10X Xenium).
    * `--min-phred-score` (int): Phred-scaled quality score cutoff (default: `20` for 10X Xenium).

    **Auxiliary Output Parameters**: 

    * `--out-transcript` (str): File name for output compressed transcript-indexed SGE file in TSV format (default: `transcripts.unsorted.tsv.gz`).
    * `--out-minmax` (str): File name for coordinate min-max values in TSV format (default: `coordinate_minmax.tsv`).
    * `--out-feature` (str): File name for compressed UMI count per gene in TSV format (default: `feature.clean.tsv.gz`).
    * `--precision-um` (int): Decimal precision for transcript coordinates; set to `0` to round to integers (default: 2).
    * `--colname-x` (str): Column name for the X-coordinate in the output SGE (default: X).
    * `--colname-y` (str): Column name for the Y-coordinate in the output SGE (default: Y).
    * `--colnames-count` (str): Comma-separated column names for expression count in the output SGE (default: count).
    * `--colname-feature-name` (str): Column name for the gene name in the output SGE(default: gene).

    **Auxiliary Environment Parameters**  
    If the binaries are already available in your system's `PATH`, you may omit these options.

    * `--gzip` (str): Path to `gzip` binary; consider `pigz -p 4` for faster processing. (Default: `gzip`)
    * `--spatula` (str): Path to `spatula` binary. (Default: `spatula`)
    * `--parquet-tools` (str): Required if `--in-parquet` is used; path to `parquet-tools` binary. (Default: `parquet-tools`)

### (Optional) Density-based Filtering

* `--filter-by-density` (flag): Enable filtering of SGE by density.
* `--out-filtered-prefix` (str): Prefix for output filtered SGE files (default: filtered).
* `--genomic-feature` (str): Genomic feature to be used for density-based filtering. Defaults to the value of `--colnames-count` if only one column name is provided.

??? note "Auxiliary Density-based Filtering Paramaters"

    * `--mu-scale` (float): Scale factor for the polygon area calculation (default: 1.0).
    * `--radius` (int): Radius for the polygon area calculation (default: 15).
    * `--quartile` (int): Quartile for the polygon area calculation (default: 2).
    * `--hex-n-move` (int): Sliding step (default: 1).
    * `--polygon-min-size` (int): The minimum polygon size (default: 500).

### (Optional) SGE Visualization
* `sge-visual` (flag): Enable SGE visualization.
* `--north-up` (flag): Enable the north-up orientation for the SGE visualization.

??? note "Auxiliary SGE Visualization Paramaters"

    * `--out-xy` (str): File name for output SGE visualization image (default: `xy.png`).
    * `--out-northup-tif` (str): File name for output north-up orientated image (default: `xy_northup.tif`).
    * `--srs` (str): If `--north-up`, define the spatial reference system (default: EPSG:3857).
    * `--resample` (str): If `--north-up`, Define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.
    * `--gdal_translate` (str): Required if `--north-up`; path to `gdal_translate` binary. (Default: `gdal_translate`)
    * `--gdalwarp` (str): Required if `--north-up`; path to `gdalwarp` binary. (Default: `gdalwarp`)

---

## Output

`cartloader` generates the following harmonized outputs:

### Unified SGE matrix

Both SGE conversion and density-based filtering generate a unified SGE matrix, consisting of:

{%
  include-markdown "../../includes/includemd_vigenettes_sgeformat.md"
%}

### SGE Images

* When `--sge-visual` is enabled, a monochrome PNG image is generated to visualize the SGE data.
* When `--north-up` is enabled, a georeferenced TIFF image is produced with a north-up orientation.