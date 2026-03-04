# Vizgen MERSCOPE Starter Tutorial

The input is an SGE dataset from the adult mouse hippocampus, extracted by masking a coronal brain section (Slice Number: 2; Replicate Number: 1; file: `detected_transcripts.csv`) from [Vizgen MERSCOPE Neuroscience Showcase](https://vizgen.com/applications/neuroscience-showcase/).

It includes steps of input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

---

## Prepare Input

### Data Access

The example data is hosted on Zenodo. Follow the commands below to download the example data.

```bash
cd $work_dir
wget  https://zenodo.org/records/17953582/files/merscope_starter.raw.tar.gz
tar --strip-components=1 -zxvf merscope_starter.raw.tar.gz
```

### File Format

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_merscope.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="merscope_hippo"                # change this to reflect your dataset name
PLATFORM="vizgen_merscope"              # platform information
SCALE=1                                 # scale from coordinate to micrometer

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

!!! info "How to Define Scaling Factors for MERSCOPE?"

    The MERSCOPE example data currently used here provides SGE in µm. Define the scaling factor from coordinate units to micrometers as 1.

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
  include-markdown "../../../includes/includemd_vigenettes_output_merscope.md"
%}
