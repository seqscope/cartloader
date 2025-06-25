# 10X Xenium Starter Tutorial

## Input Data

This tutorial uses SGE data generated with the 10x Genomics Xenium platform, and it has been cropped to a small region of the adult mouse brain for demonstration purposes.

**File Format**

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_xenium.md"
%}

**Data Access**

The example data is hosted on Zenedo ().

Follow the commands below to download the example data.

```bash
wget  
```


## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
cd $work_dir

# Unique identifier for your dataset
DATA_ID="xenium_hippo"                  # change this to reflect your dataset name
PLATFORM="10x_xenium"                   # platform information
SCALE=1                                 # coordinate to micrometer scaling factor

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

!!! info "How to Define Scaling Factors for Xenium"

    The Xenium example data currently used here provides SGE in micrometer units. Use define scaling factor from coordinate to micrometer as 1.

## SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_tsv.md"
%}

## `FICTURE` analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
%}

## `cartloader` Compilation

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
%}

## Upload to Data Repository

{%
  include-markdown "../../../includes/includemd_vigenettes_upload2aws.md"
%}
