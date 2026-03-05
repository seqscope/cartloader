# Quick Start: Run Locally

This tutorial walks you through running the `CartLoader` workflow using a minimal example dataset from the mouse hippocampus.

!!! info "Use Cases"
    This tutorial is ideal for users who want to:

    - Take full control over the environment
    - Customize workflow
    - Stay up-to-date with the latest development versions

!!! warning "Requirements"
    Users will need to:

    - Set up `CartLoader` and its dependencies locally (see [Installation](../../installation.md) guide).

---

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

---

## Prepare Input

### Data Access

The input example data is hosted on Zenodo [DOI: 10.5281/zenodo.15701393](https://doi.org/10.5281/zenodo.15701393).

Download the example data:

```bash
mkdir -p ${work_dir}/sge && cd ${work_dir}/sge

wget https://zenodo.org/records/17953582/files/seqscope_starter.std.tar.gz
tar -zxvf seqscope_starter.std.tar.gz
```

### File Format

The input is a mouse hippocampus SGE, already converted to a format compatible with `FICTURE` using [`sge_convert`](../../reference/sge_convert.md) in `CartLoader`.

{%
  include-markdown "../../../includes/includemd_vigenettes_sgeformat.md"
%}

### Define ID and Parameters

```bash
cd ${work_dir}

# Unique identifier for your dataset
DATA_ID="seqscope_hippo"                # change this to reflect your dataset name
PLATFORM="seqscope"                     # platform information

# LDA parameters
train_width=18                            # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                             # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

---

## SGE Format Conversion

The example dataset is already provided in `FICTURE`-compatible SGE format, so this conversion step is not required in this quickstart.

---

## `FICTURE` Analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
%}

---

## `CartLoader` Asset Packaging

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
%}

---

## Upload to Data Repository

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md" preserve-includer-indent=false
%}

---

## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_seqscope.md"
%}
