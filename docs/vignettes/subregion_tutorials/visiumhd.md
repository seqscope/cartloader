# 10X VisiumHD Starter Tutorial

## Input Data

The input data originates from the mouse hippocampus and is extracted from an official release of mouse brain spatial gene expression (SGE) data.

**File Format**

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_visiumhd.md"
%}


**Data Access**

The example data is hosted on Zenedo ().

Follow the commands below to download the example data.

```bash
work_dir=/path/to/work/directory
cd $work_dir
wget  https://zenodo.org/records/15701394/files/visiumhd_starter.raw.tar.gz 
tar -zxvf visiumhd_starter.raw.tar.gz  
```

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
# Unique identifier for your dataset
DATA_ID="visiumhd_hippo"                 # change this to reflect your dataset name
PLATFORM="10x_visium_hd"                 # platform information

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

!!! info "How to Define Scaling Factors for Visium HD"

    10x Visium HD includes a `scalefactors_json.json` file that provides pixel-to-micrometer scaling information. `cartloader` can directly accept this file via the `--scale-json` option and will automatically compute the appropriate scaling factor, omitting manually calculate and specify `--units-per-um`. 
    
    Alternatively, users may bypass the JSON file by directly providing a value through the `--units-per-um` option.


## SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_visiumhd.md"
%}

## `FICTURE` analysis

<!-- !! Visium HD uses --decode-scale this should be updated. -->
<!-- 
```bash
cartloader run_ficture2 \
  --makefn run_ficture2.mk \
  --main \
  --in-transcript ./sge/transcripts.unsorted.tsv.gz \
  --in-feature ./sge/feature.clean.tsv.gz \
  --in-minmax ./sge/coordinate_minmax.tsv \
  --cmap-file ${cmap} \
  --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \
  --decode-scale 2 \
  --out-dir ./ficture2 \
  --width ${train_width} \
  --n-factor ${n_factor} \
  --spatula ${spatula} \
  --ficture2 ${punkst} \
  --n-jobs 10 \
  --threads 10
``` -->

{% include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md" %}


## `cartloader` Compilation

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
%}

## Upload to Data Repository
{%
  include-markdown "../../../includes/includemd_vigenettes_upload2aws.md"
%}