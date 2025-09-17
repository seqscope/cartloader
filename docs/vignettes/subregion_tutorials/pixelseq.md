# Pixel-Seq Starter Tutorial

## Input Data

Since the Pixel-Seq publication provides SGE data only from the mouse olfactory bulb and parabrachial nucleus — neither of which includes the hippocampus — we extract a subregion from the olfactory bulb as the input for this tutorial.

**File Format**

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_pixelseq.md"
%}

**Data Access**

The example data is hosted on Zenodo.

Follow the commands below to download the example data.

```bash
work_dir=/path/to/work/directory
cd $work_dir
wget  https://zenodo.org/records/15786632/files/pixelseq_starter.raw.tar.gz 
tar --strip-components=1 -zxvf pixelseq_starter.raw.tar.gz 
```


## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
# Unique identifier for your dataset
DATA_ID="pixelseq_hippo"                # change this to reflect your dataset name
PLATFORM="pixel_seq"                    # platform information
SCALE=3.076923                        # scale from coordinate to micrometer

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

!!! info "How to Define Scaling Factors for Pixel-Seq?"

    In [Pixel-Seq publication](https://doi.org/10.1016/j.cell.2022.10.021):
    > "Because polonies have varied sizes and shapes, to maximize the feature resolution we developed a base-calling pipeline to determine the major barcode species in each pixel (0.325 * 0.325 mm2) of gel images to construct a spatial barcode map".
    
    Accordingly, we defined scale as 1/0.325 = 3.076923

## SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_tsv.md"
%}

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


## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_pixelseq.md"
%}
