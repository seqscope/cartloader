# CosMX SMI Starter Tutorial

## Input Data

The input data is from an adult mouse hippocampus, extracted by masking a coronal brain section. The original full-section 

**File Format**

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_cosmxsmi.md"
%}

**Data Access**

The example data is hosted on Zenedo ().

Follow the commands below to download the example data.

```bash
work_dir=/path/to/work/directory
cd $work_dir
wget  https://zenodo.org/records/15786632/files/cosmxsmi_starter.raw.tar.gz
tar --strip-components=1 -zxvf cosmxsmi_starter.raw.tar.gz
```

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
# Unique identifier for your dataset
DATA_ID="cosmxsmi_hippo"                # change this to reflect your dataset name
PLATFORM="cosmx_smi"                    # platform information
SCALE=$(echo 1000/120|bc -l)              # scale from coordinate to micrometer

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)

```

!!! info "How to Define Scaling Factors for CosMX SMI?"

    According to the README.html provided with the Pixel-seq dataset, each pixel has an edge length of 120 nm. To calculate the number of pixels per micrometer, use the formula: scale = 1000 / 120.

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

## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_cosmxsmi.md"
%}
