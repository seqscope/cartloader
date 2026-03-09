# CosMX SMI Starter Tutorial

This tutorial walks through a starter end-to-end workflow for CosMX SMI data using an adult mouse hippocampus subset extracted from a coronal brain section.

It includes steps of input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

---

## Prepare Input

### Data Access

The example input data is hosted on Zenodo. Follow the commands below to download it.

```bash
cd $work_dir
wget  https://zenodo.org/records/17953582/files/cosmxsmi_starter.raw.tar.gz
tar --strip-components=1 -zxvf cosmxsmi_starter.raw.tar.gz
```

### File Format

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_cosmxsmi.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="cosmxsmi_hippo"                # change this to reflect your dataset name
PLATFORM="cosmx_smi"                    # platform information
SCALE=$(echo 1000/120|bc -l)            # scale from coordinate to micrometer

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)

```

!!! question "How to define Scaling Factors for CosMX SMI?"

    According to the README.html provided with the example CosMX dataset, each pixel has an edge length of 120 nm. To calculate the number of pixels per micrometer, use the formula: scale = 1000 / 120.

---

## SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_tsv.md"
%}

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
  include-markdown "../../../includes/includemd_vigenettes_output_cosmxsmi.md"
%}
