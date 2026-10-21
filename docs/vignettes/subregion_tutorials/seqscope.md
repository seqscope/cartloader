# Seq-Scope Starter Tutorial

This tutorial uses an example SGE from *mouse hippocampus*, extracted via spatial masking from a [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0) coronal brain slice.

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step: input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Run with `run_together`

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/seqscope_starter.raw.tar.gz
tar -zxvf seqscope_starter.raw.tar.gz   # creates raw/ with barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
```

Then run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform seqscope \
      --in-dir ./raw \
      --out-dir ./out/seqscope_hippo \
      --id seqscope_hippo \
      --width 18 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

=== "Docker"

    ```bash
    docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest \
      run_together \
      --platform seqscope \
      --in-dir /data/raw \
      --out-dir /data/out/seqscope_hippo \
      --id seqscope_hippo \
      --width 18 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

The [Seq-Scope profile](../../reference/platforms/seqscope.md) reads the MEX directory given by `--in-dir` and converts its nanometer coordinates to µm. The result to deploy is `out/seqscope_hippo/cartl/seqscope_hippo/` (a `catalog.yaml` and the PMTiles it lists). Add `--dry-run` to print the commands without running them.

{%
  include-markdown "../../../includes/includemd_vigenettes_run_together_diff.md"
%}

In addition, the Seq-Scope profile counts the matrix's second count column ("GeneFull", `--icols-mtx 2`), whereas `sge_convert` below uses the first ("Gene", its default), so the counts differ slightly.

---

## Step by Step (Individual Commands)

The same analysis, one command at a time. Use it to understand each stage or to customize it.

### Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

### Input Data

Use the data downloaded [above](#run-with-run_together). Actual input formats are platform-dependent; see the other [platform tutorials](../intro.md) for theirs.

{%
  include-markdown "../../../includes/includemd_vigenettes_inputformat_seqscope.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="seqscope_hippo"                # change this to reflect your dataset name
PLATFORM="seqscope"                     # platform information
SCALE=1000                              # coordinate units per micrometer

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

!!! question "How to define Scaling Factors for Seq-Scope?"

    The latest [SeqScope](https://www.nature.com/articles/s41596-024-01065-0) with an Illumina NovaSeq 6000 uses [`NovaScope`](https://seqscope.github.io/NovaScope/) pipeline to process sequencing data. [`NovaScope`](https://seqscope.github.io/NovaScope/) defaults to generate SGE at nanometer (nm) resolution, meaning each pixel corresponds to 1 nm.

    Thus, use 1000 as scaling factor from coordinate to micrometer since 1000 nm = 1 µm.

### SGE Format Conversion

{%
  include-markdown "../../../includes/includemd_vigenettes_sge_convert_seqscope.md"
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
  include-markdown "../../../includes/includemd_vigenettes_output_seqscope.md"
%}
