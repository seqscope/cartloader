
# 🔧 Quick Start: Run Locally

This tutorial walks you through running the `CartLoader` workflow using a minimal example dataset from the mouse hippocampus.

!!! info "Use Cases"
    This tutorial is ideal for users who want to:

    - Take full control over the environment
    - Customize workflow
    - Stay up-to-date with the latest development versions 

!!! warning "Requirements"
    Users will need to:

    * Set up `CartLoader` and its dependencies locally (see [Installation](../../installation.md) guide).
    * Download the example input data (see [Input Data](#input-data))

## Input Data

The input is a mouse hippocampus SGE, already converted to a format compatible with `FICTURE` using [`sge_convert`](../../reference/sge_convert.md) in `CartLoader`.
<!-- 
!!! warning "Prepare Input SGE for `FICTURE`"

    `FICTURE` requires input in the form of a transcript-indexed SGE file in TSV format with at least: X and Y spatial coordinates, gene identifiers, and expression counts. 
    
    Because ST platforms vary widely in their data formats and metadata structures, `CartLoader` provides the [`sge_convert`](../../reference/sge_convert.md) module to convert raw SGE data into the standardized format required by `FICTURE`.
    
    For detailed platform-specific instructions on preparing compatible SGE files, see the [Vignettes](../intro.md#getting-started-per-platform). -->

**File Format**

SGE in FICTURE-compatible format includes:

{%
  include-markdown "../../../includes/includemd_vigenettes_sgeformat.md"
%}

**Data Access**

The input example data is hosted on Zenodo [DOI: 10.5281/zenodo.15701393](https://doi.org/10.5281/zenodo.15701393).

Download the example data:

```bash
work_dir=/path/to/work/directory        # path to work directory that contains the downloaded input data
cd $work_dir

wget https://zenodo.org/records/15786632/files/seqscope_starter.std.tar.gz
tar -zxvf seqscope_starter.std.tar.gz
```

---------------

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
# Unique identifier for your dataset
DATA_ID="seqscope_hippo"                # change this to reflect your dataset name
PLATFORM="seqscope"                     # platform information

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```


## `FICTURE` Analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
%}

## `CartLoader` Asset Packaging

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
%}

## Upload to Data Repository

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md" preserve-includer-indent=false
%}

---------------

## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_seqscope.md"
%}
