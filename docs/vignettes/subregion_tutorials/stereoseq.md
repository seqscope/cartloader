# Stereo-seq Starter Tutorial

The input data represents the mouse hippocampus. For demonstration purposes, we selected an adult mouse brain coronal section from the official release and extracted only the hippocampus region. The full-section data is sourced from the [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq) platform as part of the [MOSTA](https://db.cngb.org/stomics/mosta/download/) project.

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step: input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Run with `run_together`

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/stereoseq_starter.raw.tar.gz
tar --strip-components=1 -zxvf stereoseq_starter.raw.tar.gz   # extracts input.tsv.gz
```

Then run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform generic \
      --raw-transcript ./input.tsv.gz \
      --colname-transcript-x x \
      --colname-transcript-y y \
      --colname-transcript-feature geneID \
      --colname-transcript-count MIDCounts \
      --units-per-um 2 \
      --out-dir ./out/stereoseq_hippo \
      --id stereoseq_hippo \
      --width 18 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

=== "Docker"

    ```bash
    docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest \
      run_together \
      --platform generic \
      --raw-transcript /data/input.tsv.gz \
      --colname-transcript-x x \
      --colname-transcript-y y \
      --colname-transcript-feature geneID \
      --colname-transcript-count MIDCounts \
      --units-per-um 2 \
      --out-dir /data/out/stereoseq_hippo \
      --id stereoseq_hippo \
      --width 18 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

The [Stereo-seq profile](../../reference/platforms/stereoseq.md) ingests SAW `.gef` files (and needs the SAW binary), so a text GEM file like this one goes through the [generic profile](../../reference/platforms/generic.md) instead: the `--colname-transcript-*` options name its columns, and `--units-per-um 2` converts its 0.5 µm grid to µm. The result to deploy is `out/stereoseq_hippo/cartl/stereoseq_hippo/` (a `catalog.yaml` and the PMTiles it lists). Add `--dry-run` to print the commands without running them.

{%
  include-markdown "../../../includes/includemd_vigenettes_run_together_diff.md"
%}

---

## Step by Step (Individual Commands)

The same analysis, one command at a time. Use it to understand each stage or to customize it.

### Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

### Input Data

Use the data downloaded [above](#run-with-run_together).

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_stereoseq.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="stereoseq_hippo"                # change this to reflect your dataset name
PLATFORM="bgi_stereoseq"                 # platform information
SCALE=2                                  # coordinate units per micrometer

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

!!! question "How to define Scaling Factors for Stereo-seq?"

    According to the [Stereo-seq paper](https://doi.org/10.1016/j.cell.2022.04.003), the technology features spot sizes of approximately 220 nanometers in diameter with center-to-center distances of about 500 nanometers. Thus, each pixel corresponds to approximately 0.5 micrometers. The scaling factor from coordinate to µm is defined as 2.

### SGE Format Conversion

Convert the raw input to the unified SGE format. See more details in its [Reference page](../../reference/sge_convert.md). The `--csv-colname-*` options name the starter file's gene and count columns (see the note above).

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform ${PLATFORM} \
  --in-csv ./input.tsv.gz \
  --csv-colname-feature-name geneID \
  --csv-colname-count MIDCounts \
  --units-per-um ${SCALE} \
  --out-dir ./sge \
  --exclude-feature-regex '^(BLANK|Neg|Intergenic|Deprecated|Unassigned)' \
  --sge-visual \
  --n-jobs ${n_jobs}
```

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_tsv.md"
  start="<!--parameter-start-->"
  end="<!--parameter-end-->"
%}

### `FICTURE` Analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
%}

### `CartLoader` Asset Packaging

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
%}

### Upload to Data Repository

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md" preserve-includer-indent=false
%}

---

## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_stereoseq.md"
%}
