# CosMx SMI Starter Tutorial

This tutorial walks through a starter end-to-end workflow for CosMx SMI data using an adult mouse hippocampus subset extracted from a coronal brain section.

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step: input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Run with `run_together`

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/cosmxsmi_starter.raw.tar.gz
tar --strip-components=1 -zxvf cosmxsmi_starter.raw.tar.gz   # extracts input.tsv.gz
```

The [CosMx SMI profile](../../reference/platforms/cosmx_smi.md) expects a full AtoMx flat-file export (transcript, cell-metadata, and polygon files). The starter data is a single transcript table, so it goes through the [generic profile](../../reference/platforms/generic.md) with the `cosmx_smi` input preset of `sge_convert` (comma-separated; columns `x_global_px`, `y_global_px`, `target`). Write a small config file that selects the preset:

```bash
cat > cosmxsmi_starter.json <<'EOF'
{ "ingest": { "sge_platform": "cosmx_smi" } }
EOF
```

Then run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform generic \
      --raw-transcript ./input.tsv.gz \
      --config cosmxsmi_starter.json \
      --units-per-um 8.333333 \
      --out-dir ./out/cosmxsmi_hippo \
      --id cosmxsmi_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

=== "Docker"

    ```bash
    docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest \
      run_together \
      --platform generic \
      --raw-transcript /data/input.tsv.gz \
      --config /data/cosmxsmi_starter.json \
      --units-per-um 8.333333 \
      --out-dir /data/out/cosmxsmi_hippo \
      --id cosmxsmi_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

`--units-per-um 8.333333` (1000 / 120) converts the 120 nm pixels to µm. The result to deploy is `out/cosmxsmi_hippo/cartl/cosmxsmi_hippo/` (a `catalog.yaml` and the PMTiles it lists). Add `--dry-run` to print the commands without running them.

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
  include-markdown "../../../includes/includemd_vigenettes_inputformat_cosmxsmi.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="cosmxsmi_hippo"                # change this to reflect your dataset name
PLATFORM="cosmx_smi"                    # platform information
SCALE=$(echo 1000/120|bc -l)            # coordinate units per micrometer

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

!!! question "How to define Scaling Factors for CosMx SMI?"

    According to the README.html provided with the example CosMx dataset, each pixel has an edge length of 120 nm. To calculate the number of pixels per micrometer, use the formula: scale = 1000 / 120.

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
  include-markdown "../../../includes/includemd_vigenettes_output_cosmxsmi.md"
%}
