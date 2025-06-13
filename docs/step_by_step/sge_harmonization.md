# Spatial Gene Expression Matrix (SGE) Harmonization

## Overview

Spatial Digital Gene Expression (SGE) datasets vary widely in format and resolution across platforms. To enable consistent downstream analysis, `cartloader` toolkit provides a `sge_convert` module to harmonize raw SGE data by converting it into standardized transcript-based SGE in a compressed TSV format, with spatial coordinates converted to a micrometer-based unit system without altering the original resolution.

The current `sge_convert` supports standarizing SGE from sequencing-based platforms (e.g., [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0), [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq), [Pixel-seq](https://www.cell.com/cell/fulltext/S0092-8674(22)01367-8), [10x Visium HD](https://www.10xgenomics.com/platforms/visium)) and imaging-based platforms (e.g., [10x Xenium](https://www.10xgenomics.com/platforms/xenium), [Vizgen MERSCOPE](https://vizgen.com/merscope-ultra/), [CosMx SMI](https://nanostring.com/products/cosmx-spatial-molecular-imager)).

## Example Usage

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
    --icols-mtx 1 \
    --in-parquet /path/to/input/parquet/file \
    --scale-json /path/to/input/json/file \
    --exclude-feature-regex '^(BLANK|NegCon|NegPrb|mt-|Gm\d+$)' \
    --out-dir /path/to/output/dir \
    --spatula /path/to/spatul/binary \
    --sge-visual 
```

### Input SGE in TSV/CSV Format

This applied to input SGE in TSV/CSV format from platforms including 10X Xenium, StereoSeq, Cosmx SMI, MERSCOPE, Pixel-seq. Below is an example converting SGE from StereoSeq.

```bash
cartloader sge_convert \
    --platform bgi_stereoseq \
    --in-csv /path/to/input/csv/file \
    --units-per-um 28.75 \
    --out-dir /path/to/output/dir \
    --colnames-count count \
    --exclude-feature-regex '^(BLANK|NegCon|NegPrb|mt-|Gm\d+$)' \
    --filter-by-density \
    --out-filtered-prefix filtered \
    --genomic-feature count \
    --sge-visual \
    --north-up \
    --spatula /path/to/spatul/binary 
```

---
## Actions

### SGE Conversion

See more parameters using `cartloader sge_convert --help`.

#### Required Parameters

* `--platform`: Source platform to infer the input file format.
* `--in-mex`: Path to the input SGE directory in MEX format. Required for MEX-formatted data (e.g., 10X Visium HD, SeqScope).
* `--in-csv`: Path to the input SGE file in CSV or TSV format. Required for CSV/TSV-formatted data (e.g., 10X Xenium, StereoSeq, CosMX SMI, MERSCOPE, Pixel-seq, NovaST).
* `--in-parquet`: Path to the input Parquet file containing spatial coordinates, if available. Typically named `tissue_positions.parquet` (10X Visium HD).
* `--scale-json`: Path to the scale JSON file used to compute `--units-per-um`, if available. Typically named `scalefactors_json.json` (10X Visium HD).
* `--units-per-um`: Coordinate units per micrometer (default: 1.00). Skip if `--scale-json` is provided.
* `--out-dir`: Path to the output directory.

#### Auxiliary Input MEX Parameters

* `--icols-mtx`: Comma-separated, 1-based indices of the target genomic features among the count columns in the input matrix file. (Default: 1)
* `--colnames-count`: Comma-separated output column names for the specified genomic features. (Default: count). The number of names specified by `--colnames-count` must match the number of indices provided in `--icols-mtx`.

#### Auxiliary Input CSV/TSV Parameters

* `--csv-comment`: If enabled, lines starts with `#` will be skipped (default: `False` for 10X Xenium, StereoSeq, CosMx SMI, MERSCOPE, and Pixel-seq; `True` for NovaST).
* `--csv-delim`: Delimiter for the input file (default: `","` for 10X Xenium, CosMx SMI, and MERSCOPE; `"\t"` for StereoSeq, Pixel-seq, and NovaST).
* `--csv-colname-x`: Column name for X coordinates (default: `x_location` for 10X Xenium; `x` for StereoSeq and NovaST; `x_local_px` for CosMx SMI; `global_x` for MERSCOPE; `xcoord` for Pixel-seq).
* `--csv-colname-y`: Column name for Y coordinates (default: `y_location` for 10X Xenium; `y` for StereoSeq and NovaST; `y_local_px` for CosMx SMI; `global_y` for MERSCOPE; `ycoord` for Pixel-seq).
* `--csv-colnames-count`: Comma-separated column names for expression count. If not provided, a count of 1 per transcript (default: `MIDCounts` for StereoSeq; `MIDCount` for NovaST).
* `--csv-colname-feature-name`: Column name for gene name (default: `feature_name` for 10X Xenium; `geneID` for StereoSeq; `target` for CosMx SMI; `gene` for MERSCOPE; `geneName` for Pixel-seq; `geneID` for NovaST).
* `--csv-colname-feature-id`: Column name for gene id.
* `--csv-colnames-others`: Columns names to keep.
* `--csv-colname-phredscore`: Column name for Phred-scaled quality value estimating the probability of incorrect calls (default: `qv` for 10X Xenium).
* `--min-phred-score`: Phred-scaled quality score cutoff (default: `20` for 10X Xenium).

#### Auxiliary Output Parameters

* `--out-transcript`: File name for output compressed transcript-indexed SGE file in TSV format (default: `transcripts.unsorted.tsv.gz`).
* `--out-minmax`: File name for coordinate min-max values in TSV format (default: `coordinate_minmax.tsv`).
* `--out-feature`: File name for compressed UMI count per gene in TSV format (default: `feature.clean.tsv.gz`).
* `--precision-um`: Decimal precision for transcript coordinates; set to `0` to round to integers (default: 2).
* `--colname-x`: Column name for the X-coordinate in the output SGE (default: X).
* `--colname-y`: Column name for the Y-coordinate in the output SGE (default: Y).
* `--colnames-count`: Comma-separated column names for expression count in the output SGE (default: count).
* `--colname-feature-name`: Column name for the gene name in the output SGE(default: gene).
* `--colname-feature-id`: Column name for the gene ID in the output SGE; required only when `--csv-colname-feature-id` is defined (default: None).

### Optional Density-based Filtering

* `--filter-by-density`: Enable filtering of SGE by density.
* `--out-filtered-prefix`: Prefix for output filtered SGE files (default: filtered).
* `--genomic-feature`: Genomic feature to be used for density-based filtering. Defaults to the value of `--colnames-count` if only one column name is provided.
* `--mu-scale`: Scale factor for the polygon area calculation (default: 1.0).
* `--radius`: Radius for the polygon area calculation (default: 15).
* `--quartile`: Quartile for the polygon area calculation (default: 2).
* `--hex-n-move`: Sliding step (default: 1).
* `--polygon-min-size`: The minimum polygon size (default: 500).

### Optional Gene Filtering

* `--include-feature-list` and `--exclude-feature-list`: A file containing a list of input genes to be included or excluded.
* `--include-feature-regex` and `--exclude-feature-regex`: A regex pattern of feature/gene names to be included or excluded.
* `--include-feature-type-regex`: A regex pattern of feature/gene type to be included.
* `--csv-colname-feature-type`: If `--include-feature-type-regex`, column name in the input TSV/CSV that contains gene type information.
* `--feature-type-ref`, `--feature-type-ref-delim`, `--feature-type-ref-colidx-name`, and `--feature-type-ref-colidx-type`: If `--include-feature-type-regex`, define the file path, delimiter, and column indices for gene names and types of an additional reference file to provide gene type information.

### Optional SGE Visualization

* `sge-visual`: Enable SGE visualization.
* `--out-xy`: Output for SGE visualization image (default: xy.png)
* `--north-up`: Enable the north-up orientation for the SGE visualization.
* `--out-northup-tif`:  Output for SGE visualization after north-up orientation. (default: `xy_northup.tif`).
* `--srs`: If `--north-up`, define the spatial reference system (default: EPSG:3857).
* `--resample`: If `--north-up`, Define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.

---

## Output

Cartloader generates the following harmonized outputs:

### Unified SGE matrix

Both SGE conversion and density-based filtering generate a unified spatial gene expression (SGE) matrix, consisting of:

* A compressed transcript-indexed SGE file in TSV format.
* A TSV file for min and max X and Y coordinates.
* A TSV file collects UMI counts on a per-gene basis.

### SGE Images

* When `--sge-visual` is enabled, a monochrome PNG image is generated to visualize the SGE data.
* When `--north-up` is enabled, a georeferenced TIFF image is produced with a north-up orientation.