# StereoSeq Starter Tutorial

## Input Data

The input data represents mouse hippocampus. For demonstration purposes, we selected the adult mouse brain coronal section from the official release and extracted only the hippocampus region. The full-section data sourced from the [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq) platform, as part of the [MOSTA](https://db.cngb.org/stomics/mosta/download/) project.

**File Format**

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_stereoseq.md"
%}

**Data Access**

The example data is hosted on Zenedo ().

Follow the commands below to download the example data.

```bash
work_dir=/path/to/work/directory
cd $work_dir
wget  https://zenodo.org/records/15786632/files/stereoseq_starter.raw.tar.gz 
tar --strip-components=1 -zxvf stereoseq_starter.raw.tar.gz  
```

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
# Unique identifier for your dataset
DATA_ID="stereoseq_hippo"                # change this to reflect your dataset name
PLATFORM="bgi_stereoseq"                 # platform information
SCALE=2                                  # scale from coordinate to micrometer

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

!!! info "How to Define Scaling Factors for StereoSeq?"

    According to the [StereoSeq paper](https://doi.org/10.1016/j.cell.2022.04.003), the technology features spot sizes of approximately 220 nanometers in diameter with center-to-center distances of about 500 nanometers. Thus, each pixel corresponds to approximately 0.5 micrometers. The scaling factor from coordinate to um is defined as 2.

## SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_tsv.md"
%}

## `FICTURE` Analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
%}

## `cartloader` Compilation

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
%}

## Upload to Data Repository
{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md"
%}

## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_stereoseq.md"
%}
