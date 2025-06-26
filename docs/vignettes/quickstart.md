# A Quick Start

This page provides a quick-start example with a small dataset from mouse hippocampus to help users verify their installation and get familiar with `cartloader` and the `cartostore` workflow.

The workflow consists of three major steps:

1.	`FICTURE` Analysis — computes factor embeddings and visualizations via `punkst`.
2.	`cartloader` Compilation — generates web-compatible tiles.
3.	`AWS` Upload — places results into an AWS S3 bucket for sharing via `cartostore`.


!!! warning "Requirement"

    Before getting started, please ensure that cartloader and all prerequisites are installed (see [Installation](../installation.md) guide).

## Input Data

This tutorial used an SGE representing the mouse hippocampus as input data. The input SGE file was prepared in a format compatible with `FICTURE` using the [`sge_convert`](../reference/sge_convert.md) function from the `cartloader` package.

!!! warning "Prepare Input SGE for `FICTURE`"

    `FICTURE` requires input in the form of a transcript-indexed SGE file in TSV format with at least: X and Y spatial coordinates, gene identifiers, and expression counts. 
    
    Because spatial transcriptomics (ST) platforms vary widely in their data formats and metadata structures, the `cartloader` package provides the [`sge_convert`](../reference/sge_convert.md) function to convert raw SGE data into the standardized format required by `FICTURE`.
    
    For detailed platform-specific instructions on preparing compatible SGE files, see the [Vignettes](../vignettes/intro.md#getting-started-per-platform).

**File Format**

The example SGE includes the following files:

??? "`transcripts.unsorted.tsv.gz`: transcript-indexed SGE in TSV"

    ```text
    X        Y        gene     count
    295.29   1422.35  Myo3a    0
    1745.54  1110.72  Med14    1
    1745.54  1110.72  Ntpcr    1
    ```

    * `X`: X coordinates in um
    * `Y`: Y coordinates in um
    * `gene`: gene symbols
    * `count`: expression count per pixel per gene

??? "`feature.clean.tsv.gz`: (optional) UMI counts on a per-gene basis in TSV"

    ```text
    gene           gene_id             count
    Gm29155        ENSMUSG00000100764  1
    Pcmtd1         ENSMUSG00000051285  431
    Gm26901        ENSMUSG00000097797  1
    ```
    * `gene`: gene symbols
    * `gene_id`: gene IDs
    * `count`: expression count per gene

??? "`coordinate_minmax.tsv`: (optional) X Y min/max coordinates"

    ```text
    xmin	0.14
    xmax	2359.90
    ymin	0.23
    ymax	1439.95
    ```

    * `xmin` `xmax`: min and max X coordinates in um
    * `ymin` `ymax`: min and max Y coordinates in um


**Data Access**

The example data is hosted on Zenedo (10.5281/zenodo.15701394).

Follow the commands below to download the example data.

```bash
work_dir=/path/to/work/directory
cd $work_dir

wget  https://zenodo.org/records/15701394/files/seqscope_starter.std.tar.gz
tar -zxvf seqscope_starter.std.tar.gz
```

## Set Up the Environment

{%
  include-markdown "../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
cd $work_dir

# Unique identifier for your dataset
DATA_ID="seqscope_hippo"                # change this to reflect your dataset name
PLATFORM="seqscope"                     # platform information

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

## `FICTURE` analysis

{%
  include-markdown "../../includes/includemd_vigenettes_run_ficture2.md"
%}

## `cartloader` Compilation

{%
  include-markdown "../../includes/includemd_vigenettes_run_cartload2.md"
%}

## Upload to Data Repository

{%
  include-markdown "../../includes/includemd_vigenettes_upload2aws.md"
%}
