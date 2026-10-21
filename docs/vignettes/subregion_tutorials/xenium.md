# 10x Xenium Starter Tutorial

This tutorial uses SGE data generated with the 10x Genomics Xenium platform, and it has been cropped to a small region of the adult mouse brain for demonstration purposes.

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step: input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Run with `run_together`

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/xenium_starter.raw.tar.gz
tar --strip-components=1 -zxvf xenium_starter.raw.tar.gz   # extracts input.tsv.gz
```

The [10x Xenium profile](../../reference/platforms/xenium.md) carries the `cell_id`, `z_location`, `overlaps_nucleus`, and `nucleus_distance` columns through ingest. This starter file predates `nucleus_distance`, so write a small config file that leaves it out:

```bash
cat > xenium_starter.json <<'EOF'
{ "ingest": { "csv_colnames_others": ["cell_id", "z_location", "overlaps_nucleus"] } }
EOF
```

Then run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform 10x_xenium \
      --raw-transcript ./input.tsv.gz \
      --config xenium_starter.json \
      --out-dir ./out/xenium_hippo \
      --id xenium_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

=== "Docker"

    ```bash
    docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest \
      run_together \
      --platform 10x_xenium \
      --raw-transcript /data/input.tsv.gz \
      --config /data/xenium_starter.json \
      --out-dir /data/out/xenium_hippo \
      --id xenium_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

`--raw-transcript` names the transcript file directly, since the starter data is a single file rather than a Xenium Ranger output directory (`--in-dir`). The result to deploy is `out/xenium_hippo/cartl/xenium_hippo/` (a `catalog.yaml` and the PMTiles it lists). Add `--dry-run` to print the commands without running them.

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
  include-markdown "../../../includes/includemd_vigenettes_inputformat_xenium.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="xenium_hippo"                  # change this to reflect your dataset name
PLATFORM="10x_xenium"                   # platform information
SCALE=1                                 # coordinate units per micrometer

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

!!! question "How to define Scaling Factors for Xenium?"

    The Xenium example data currently used here provides SGE in micrometer units. Define the scaling factor from coordinate units to micrometers as 1.

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
  include-markdown "../../../includes/includemd_vigenettes_output_xenium.md"
%}
