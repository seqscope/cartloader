# Multi‑Sample Batch Analysis Tutorial

## Overview

This tutorial demonstrates how to run multi‑sample spatial factor analysis and packaging using `run_ficture2_multi` and `run_cartload2_multi`. You’ll prepare a sample list, generate Makefiles, execute them with `make`, and optionally dry‑run and upload results.

Use this when training shared spatial factors across multiple samples and generating per‑sample outputs in a single, parallelizable workflow.

---
## Input Data

**File Format**

TODO:

**Data Access**

Download the source data from Zenodo.

Follow the commands below to download the example data.

```bash
work_dir=/path/to/work/directory
cd $work_dir

wget https://zenodo.org/records/15127709/files/FB080_O1a.zip?download=1
unzip FB080_O1a.zip

wget https://zenodo.org/records/15127709/files/FB080_O1b.zip?download=1
unzip FB080_O1b.zip

wget https://zenodo.org/records/15127709/files/FB080_O1c.zip?download=1
unzip FB080_O1c.zip

wget https://zenodo.org/records/15127709/files/FB080_O1d.zip?download=1
unzip FB080_O1d.zip
```

---
## Set Up the Environment

{% 
  include-markdown "../../includes/includemd_vigenettes_setupenv.md"
%}


Define data ID and analysis parameters:

```bash
# Unique identifier for your collection
COLLECTION_ID="walsh2025-human-cortex-fb080-O1" # change this to reflect your dataset name
PLATFORM="generic"                      # platform information
SCALE=1                                 # coordinate to micrometer scaling factor

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=48,96                           # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

Prepare a working directory:

```bash
work_dir=/path/to/work/dir/${COLLECTION_ID}
mkdir -p ${work_dir}
cd ${work_dir}
```

---

## SGE Format Conversion

For each sample, convert its transcript‑indexed SGE file to be the `CartLoader` format.

TODO - add the example code.

## Prepare Input List

Create a tab-separataed TSV file (e.g., `input.tsv`) with one sample per line to provide the information for each sample. The TSV file should have no header and with two-to-four columns.

Example Input List: `input.tsv`
```text
FB080_O1a	/path/to/FB080_O1a/transcripts.tsv.gz
FB080_O1b	/path/to/FB080_O1b/transcripts.tsv.gz
FB080_O1c	/path/to/FB080_O1c/transcripts.tsv.gz
FB080_O1d	/path/to/FB080_O1d/transcripts.tsv.gz
```

* `1st Column` (required, str): Unique dataset identifier (avoid whitespace; prefer `-` than `_`)
* `2nd Column` (required, str): Path to the transcript‑indexed SGE file in `CartLoader` format, generated from [SGE conversion](../../reference/sge_convert.md)).
* `3rd Column`  (optional, str): Optional dataset title. Quote the value when whitespace is involved.
* `4th Column`  (optional, str): Optional dataset description. Quote the value when whitespace is involved.

## Multi‑Sample FICTURE Analysis

Trains LDA models cross samples per LDA parameter setting (width and n-factor), decode per sample, and summarize per‑sample JSON manifests summarizing results.

Outputs are written under `ficture2/<sample_id>/`. See details in [run_ficture2_multi](../../reference/run_ficture2_multi.md).

```bash
cartloader run_ficture2_multi \
  --in-list ./input.tsv \
  --out-dir ./ficture2 \
  --width 18 \
  --n-factor 48,96 \
  --exclude-feature-regex "^(Blank-.*$)" \
  --redo-merge-units \
  --ficture2 ${PUNKST} \
  --spatula ${SPATULA} \
  --threads 12 \
  --n-jobs 10 

```

---
## Multi‑Sample Asset Packaging

Package per‑sample outputs from `ficture2/<sample_id>/` into web‑ready PMTiles and catalogs.

Outputs are written under `cartload2/<sample_id>/`. See details in [run_cartload2_multi](../../reference/run_cartload2_multi.md).

```bash
cartloader run_cartload2_multi \
  --in-list ./input.tsv \
  --fic-dir ./ficture2 \
  --out-dir ./cartload2 \
  --spatula ${SPATULA} \
  --tippecanoe ${TIPPECANOE} \
  --threads 12 \
  --n-jobs 10 
```

---

## Upload to AWS

You can upload the generated CartLoader outputs to your designated AWS S3 directory by a single sample or the whole collection. For full details, see [upload_aws](../../reference/upload_aws.md).

=== "Upload a collection"

    ```bash
    AWS_DIR=s3://your-bucket/${COLLECTION_ID} # Recommend to use COLLECTION_ID as directory name, this will create a subdirectory for each sample

    cartloader upload_aws \
      --in-dir ./cartload2 \
      --s3-dir ${AWS_DIR} \
      --in-list ./input.tsv 
    ```

=== "Upload a single sample"

    ```bash
    AWS_DIR=s3://your-bucket/${COLLECTION_ID} # Recommend to use COLLECTION_ID as directory name, this will create a subdirectory for each sample

    cartloader upload_aws \
      --in-dir ./cartload2/FB080_O1a \
      --s3-dir ${AWS_DIR}/fb080-01a
    ```

---
## Output Summary

See reference details: [run_ficture2_multi](../../reference/run_ficture2_multi.md#output) and[run_cartload2_multi](../../reference/run_cartload2_multi.md#output).

### Spatial Factor Inference Per Sample

Below is an example of spatial factor inference results from `FICTURE` using a training width of 18, 12 factors, a fit width of 18, and an anchor resolution of 6.

![FICTURE](../images/multisample_vigenettes/FB080_O1a.t18_f96_p18_a6.png)
![cmap](../images/multi-sample.t18_f96.cmap.png)

{{ read_csv('../tabs/FB080_O1a.t18_f96_p18_a6.factor.info.tsv',sep = '\t') }}

### Packed SGE and Spatial Factor Outputs Per Sample

TO-DO: path to those assets will be provided to serve the output