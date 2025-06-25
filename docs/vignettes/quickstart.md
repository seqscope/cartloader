# A Quick Start

This page provides a quick-start example with a small dataset to help users verify their installation and get familiar with `cartloader` and the `cartostore` workflow.

The workflow consists of four major steps:

1.	SGE Format Conversion — prepares raw SGE into a harmonized format.
2.	FICTURE Analysis — computes factor embeddings and visualizations via `punkst`.
3.	Cartloader Compilation — generates web-compatible tiles using `tippecanoe` and `pmtiles`.
4.	AWS Upload — places results into an AWS S3 bucket for sharing via `cartostore`.


!!! warning "Requirement"

    Before getting started, please ensure that cartloader and all prerequisites are installed by following the [Installation](../installation.md) guide.

## Input Data

!!! info "Input Data Requirement:"
    Input data should be a transcript-indexed SGE containing at least:

    * Spatial coordinates (X coordinates, Y coordinates)
    * Feature metadata (gene symbols)
    * Expression Counts

**Example Dataset**

This tutorial uses an example SGE from *mouse hippocampus*, extracted via spatial masking from a [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0) coronal brain slice.

**File Format**

Actual input formats are platform-dependent. Please refer to the [Vignettes](../vignettes/intro.md) for detailed input specifications by each platform.

{%
  include-markdown "../../includes/includemd_vigenettes_inputformat_seqscope.md"
%}

**Data Access**

The example data is hosted on Zenedo (10.5281/zenodo.15701394).

Follow the commands below to download the example data.

```bash
wget  https://zenodo.org/api/records/15701394/draft/files/seqscope_mousebrain_subset.raw.tar.gz/content 
tar zxvf seqscope_mousebrain_subset.raw.tar.gz
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
SCALE=1000                            # scale from coordinate to micrometer

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

!!! info "How to Define Scaling Factors for Seq-Scope"

    The latest [SeqScope](https://www.nature.com/articles/s41596-024-01065-0) with an Illumina NovaSeq 6000 uses [`NovaScope`](https://seqscope.github.io/NovaScope/) pipeline to process sequencing data. [`NovaScope`](https://seqscope.github.io/NovaScope/) defaults to generate SGE at nanometer (nm) resolution, meaning each pixel corresponds to 1 nm.

    Thus, use 1000 as scaling factor from coordinate to micrometer since 1000 nm = 1 µm.


## SGE Format Conversion

{%
  include-markdown "../../includes/includemd_vigenettes_sge_convert_seqscope.md"
%}

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
