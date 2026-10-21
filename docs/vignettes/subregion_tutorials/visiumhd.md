# 10x Visium HD Starter Tutorial

The input data originates from the mouse hippocampus and is extracted from an official release of mouse brain SGE data.

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step: input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Run with `run_together`

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/visiumhd_starter.raw.tar.gz
tar -zxvf visiumhd_starter.raw.tar.gz   # creates raw/ with the 2 µm MEX files, tissue_positions.parquet, scalefactors_json.json
```

The [10x Visium HD profile](../../reference/platforms/visium_hd.md) expects a Space Ranger `outs/` directory, where these files sit under `binned_outputs/square_002um/`. The starter data keeps them flat in `raw/`, so write a small config file that points the profile at them (paths are relative to `--in-dir`):

```bash
cat > visiumhd_starter.json <<'EOF'
{
  "ingest": {
    "inputs": {
      "--in-mex": ".",
      "--pos-parquet": "tissue_positions.parquet",
      "--scale-json": "scalefactors_json.json"
    }
  }
}
EOF
```

Then run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform 10x_visium_hd \
      --in-dir ./raw \
      --config visiumhd_starter.json \
      --out-dir ./out/visiumhd_hippo \
      --id visiumhd_hippo \
      --width 18 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

=== "Docker"

    ```bash
    docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest \
      run_together \
      --platform 10x_visium_hd \
      --in-dir /data/raw \
      --config /data/visiumhd_starter.json \
      --out-dir /data/out/visiumhd_hippo \
      --id visiumhd_hippo \
      --width 18 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

The result to deploy is `out/visiumhd_hippo/cartl/visiumhd_hippo/` (a `catalog.yaml` and the PMTiles it lists). Add `--dry-run` to print the commands without running them. With a full Space Ranger `outs/` directory, drop `--config`: the profile then also imports the segmented cells, the 8 µm and 16 µm bins, and the H&E image when present.

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
  include-markdown "../../../includes/includemd_vigenettes_inputformat_visiumhd.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="visiumhd_hippo"                 # change this to reflect your dataset name
PLATFORM="10x_visium_hd"                 # platform information

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

!!! question "How to define Scaling Factors for Visium HD?"

    10x Visium HD includes a `scalefactors_json.json` file that provides pixel-to-micrometer scaling information. `CartLoader` can directly accept this file via the `--scale-json` option and automatically compute the appropriate scaling factor, so you do not need to manually calculate and specify `--units-per-um`.
    
    Alternatively, users may bypass the JSON file by directly providing a value through the `--units-per-um` option.

### SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_visiumhd.md"
%}

### `FICTURE` Analysis

Compute spatial factors using `punkst` (FICTURE2). See more details on the [Reference page](../../reference/run_ficture2.md). For Visium HD, `--decode-scale 2` draws the factor image at 2 µm per pixel, the bin size.

```bash
cartloader run_ficture2 \
  --makefn run_ficture2.mk \
  --main \
  --in-transcript ./sge/transcripts.unsorted.tsv.gz \
  --in-feature ./sge/feature.clean.tsv.gz \
  --in-minmax ./sge/coordinate_minmax.tsv \
  --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \
  --decode-scale 2 \
  --out-dir ./ficture2 \
  --width ${train_width} \
  --n-factor ${n_factor} \
  --n-jobs ${n_jobs} \
  --threads ${n_jobs}
```

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
  start="<!--parameter-start-->"
  end="<!--parameter-end-->"
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
  include-markdown "../../../includes/includemd_vigenettes_output_visiumhd.md"
%}
