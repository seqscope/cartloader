# Human Cortex Multi‑Sample Spatial Factor Analysis Tutorial

## Overview

This tutorial demonstrates how to run multi‑sample FICTURE analysis and package results using `run_ficture2_multi` and `run_cartload2_multi`.

!!! info "Why use Multi‑Sample FICTURE Analysis?"

    **`CartLoader` supports analyzing ≥2 samples in two ways**:

    - Multi‑sample FICTURE analysis ([`run_ficture2_multi`](../reference/run_ficture2_multi.md)): Jointly learns spatial factors across all samples and writes per‑sample outputs in one parallelizable run.
    - SGE stitch + single‑sample analysis ([`sge_stitch`](../reference/sge_addon.md) → [`run_ficture2`](../reference/run_ficture2.md)): Stitch multiple SGEs into a single mosaic, then train one model on that mosaic.

    **What to expect**

    - Shared factors/comparability: `run_ficture2_multi` learns a cohort‑wide latent basis and returns per‑sample decodes for direct comparison. The stitch approach yields a single model over the merged mosaic; useful when you need a unified coordinate system (e.g., tiling adjacent sections).
    - Efficiency and scale: `run_ficture2_multi` fits once for the cohort and decodes per sample, avoiding repeated runs and post‑hoc alignment. Stitching can be simpler for mosaics but often increases I/O and memory due to very large merged files.

    **Recommendation:**

      - Prefer `run_ficture2_multi` for most cohorts for clean per‑sample outputs and better computational efficiency; use stitching when a single shared coordinate frame is required.
    - If you choose stitching, plan for higher resource usage (RAM, disk, and I/O). Large mosaics can be slow to generate and train on, and may require substantially more memory and temporary storage than per‑sample runs.

---
## Input Data

This tutorial uses a series of four human cortex ST datasets from [Walsh et al. Nature 2025](https://www.nature.com/articles/s41586-025-09010-1), generated using MERFISH.

**File Format**

!!! info "`detected_transcripts.csv.gz`"

      ```
      ,barcode_id,global_x,global_y,global_z,x,y,fov,gene,transcript_id,cell_id
      36,0,154.75854,8197.335,0.0,1790.1895,303.81055,0,HS3ST1,ENST00000002596,224295401131100634
      53,0,157.40007,8211.97,0.0,1814.648,439.31885,0,HS3ST1,ENST00000002596,-1
      76,0,154.28473,8230.835,0.0,1785.8022,614.0,0,HS3ST1,ENST00000002596,224295401131100700
      ```

      * Column 1: Unique numeric index for each transcript within a field of view (non-consecutive, ascending).
      * `barcode_id`: Zero-based index of the transcript barcode in the codebook; forms a composite key with `fov`.
      * `global_x`: Transcript x coordinates (µm) in the experimental region; may be negative due to alignment.
      * `global_y`: Transcript y coordinates (µm) in the experimental region; may be negative due to alignment.
      * `global_z`: Zero‑based z‑position index.
      * `x`: The x-coordinate of the transcript (µm), within the coordinate space of the field of view.
      * `y`: The y-coordinate of the transcript (µm), within the coordinate space of the field of view.
      * `fov`: Zero-based field of view index; forms a composite key with `barcode_id`.
      * `gene`: Gene name.
      * `transcript_id`: Unique identifier for the transcript.
      * `cell_id`: Unique identifier for cell.

**Data Access**

Follow the commands below to download the source data.

```bash
work_dir=/path/to/work/directory
mkdir -p ${work_dir}/raw
cd ${work_dir}/raw

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
cd ${work_dir}
```

---

## SGE Format Conversion

For each sample, convert its transcript‑indexed SGE file to be the `CartLoader` format.

Below use FB080_O1a as an example.

```bash
sample_id=FB080_O1a
mkdir -p ./sge/${sample_id}

cartloader sge_convert \
  --in-csv ./raw/${sample_id}/detected_transcripts.csv.gz \
  --platform generic \
  --out-dir ./sge/${sample_id} \
  --csv-delim "," \
  --csv-colname-x global_x \
  --csv-colname-y global_y \
  --csv-colname-feature-name gene \
  --csv-colnames-others cell_id,global_z \
  --sge-visual 
```

Alternatively, if you don’t need SGE visualization or density/feature‑based filtering, you can directly generate a CartLoader‑compatible SGE with fixed column names (the raw SGE is already in micrometers):

```bash
sample_id=FB080_O1a
mkdir -p ./sge/${sample_id}

(echo -e "X\tY\tgene\tcount\tcell_id\tZ";gzip -cd ./raw/${sample_id}/detected_transcripts.csv.gz | tail -n +2| tr , '\t' | perl -lane 'print join("\t",$F[2],$F[3],$F[8],1,$F[10],$F[4]);';) | gzip > ./sge/${sample_id}/transcripts.unsorted.tsv.gz
```

## Prepare Input List

Create a tab‑separated TSV file (e.g., `input.tsv`) with one sample per line. The TSV should have no header and contain two to four columns.

!!! info "Example Input List: `input.tsv`"
      ```text
      FB080_O1a	/path/to/FB080_O1a/transcripts.unsorted.tsv.gz
      FB080_O1b	/path/to/FB080_O1b/transcripts.unsorted.tsv.gz
      FB080_O1c	/path/to/FB080_O1c/transcripts.unsorted.tsv.gz
      FB080_O1d	/path/to/FB080_O1d/transcripts.unsorted.tsv.gz
      ```

      * `1st Column` (required, str): Unique dataset identifier (avoid whitespace; prefer `-` to `_`).
      * `2nd Column` (required, str): Path to the transcript‑indexed SGE file in `CartLoader` format, generated from [SGE Format Conversion](#sge-format-conversion).
      * `3rd Column`  (optional, str): Dataset title. Quote the value when whitespace is involved.
      * `4th Column`  (optional, str): Dataset description. Quote the value when whitespace is involved.

## Multi‑Sample FICTURE Analysis

Trains LDA models across samples for each parameter setting (width and n‑factor), decodes per sample, and writes per‑sample JSON manifests summarizing results.

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

Upload the generated `CartLoader` outputs from above to your designated AWS S3 directory — either the whole collection at once or each sample individually. For full details, see [upload_aws](../../reference/upload_aws.md).

=== "Upload a collection"

    ```bash
    AWS_DIR=s3://your-bucket/${COLLECTION_ID} # Recommend to use COLLECTION_ID as directory name, this will create a subdirectory for each sample

    cartloader upload_aws \
      --in-dir ./cartload2 \
      --s3-dir ${AWS_DIR} \
      --in-list ./input.tsv 
    ```

=== "Upload a single sample"

    Below use FB080_O1a as an example.

    ```bash
    AWS_DIR=s3://your-bucket/${COLLECTION_ID} # Recommend to use COLLECTION_ID as directory name, this will create a subdirectory for each sample

    cartloader upload_aws \
      --in-dir ./cartload2/FB080_O1a \
      --s3-dir ${AWS_DIR}/fb080-01a
    ```

---
## Output Summary

<!-- See reference details: [run_ficture2_multi](../../reference/run_ficture2_multi.md#output) and [run_cartload2_multi](../../reference/run_cartload2_multi.md#output).

### Spatial Factor Inference Per Sample

Below is an example of spatial factor inference results from `FICTURE`, trained at width 18 with 96 factors, decoded at width 18, and anchor resolution 6.

![FICTURE](../images/multisample_vigenettes/FB080_O1a.t18_f96_p18_a6.png)
![cmap](../images/multi-sample.t18_f96.cmap.png)

{{ read_csv('../tabs/FB080_O1a.t18_f96_p18_a6.factor.info.tsv',sep = '\t') }}

### Packed SGE and Spatial Factor Outputs Per Sample

TODO: Provide paths to these assets to serve the outputs. -->


<div class="grid cards generic" markdown>

-   ![Interactive Exploration](placeholder)

    ---

    #### View/Explore

    The output are available in both CartoScope. 

    [Explore in CartoScope](http://localhost:5173/datasets?collections=human+cerebral+cortex+development%2C+gestational+week+20+%E2%80%93+occipital+cortex){ .md-button .md-button--primary .button-tight-small }

</div>

See more details of output at the Reference pages for [run_ficture2](../docs/reference/run_ficture2.md) and [run_cartload2](../docs/reference/run_cartload2.md).

