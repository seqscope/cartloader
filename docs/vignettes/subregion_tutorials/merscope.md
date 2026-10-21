# Vizgen MERSCOPE Starter Tutorial

The input is an SGE dataset from the adult mouse hippocampus, extracted by masking a coronal brain section (Slice Number: 2; Replicate Number: 1; file: `detected_transcripts.csv`) from [Vizgen MERSCOPE Neuroscience Showcase](https://vizgen.com/applications/neuroscience-showcase/).

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step: input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Run with `run_together`

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/merscope_starter.raw.tar.gz
tar --strip-components=1 -zxvf merscope_starter.raw.tar.gz   # extracts input.tsv.gz
```

Then run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform merfish \
      --raw-transcript ./input.tsv.gz \
      --out-dir ./out/merscope_hippo \
      --id merscope_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

=== "Docker"

    ```bash
    docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest \
      run_together \
      --platform merfish \
      --raw-transcript /data/input.tsv.gz \
      --out-dir /data/out/merscope_hippo \
      --id merscope_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

The [MERSCOPE profile](../../reference/platforms/merscope.md) (`--platform merfish`) reads the transcript file named by `--raw-transcript`. The result to deploy is `out/merscope_hippo/cartl/merscope_hippo/` (a `catalog.yaml` and the PMTiles it lists). Add `--dry-run` to print the commands without running them.

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
  include-markdown "../../../includes/includemd_vigenettes_inputformat_merscope.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="merscope_hippo"                # change this to reflect your dataset name
PLATFORM="vizgen_merscope"              # platform information
SCALE=1                                 # coordinate units per micrometer

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

!!! question "How to define Scaling Factors for MERSCOPE?"

    The MERSCOPE example data currently used here provides SGE in µm. Define the scaling factor from coordinate units to micrometers as 1.

### SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_tsv.md"
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
  include-markdown "../../../includes/includemd_vigenettes_output_merscope.md"
%}
