# Pixel-Seq Starter Tutorial

Since the Pixel-Seq publication provides SGE data only from the mouse olfactory bulb and parabrachial nucleus — neither of which includes the hippocampus — we extract a subregion from the olfactory bulb as the input for this tutorial.

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step: input preparation, SGE format conversion, FICTURE analysis, asset packaging, and data upload.

---

## Run with `run_together`

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/pixelseq_starter.raw.tar.gz
tar --strip-components=1 -zxvf pixelseq_starter.raw.tar.gz   # extracts input.tsv.gz
```

Then run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform generic \
      --raw-transcript ./input.tsv.gz \
      --colname-transcript-x xcoord \
      --colname-transcript-y ycoord \
      --colname-transcript-feature geneName \
      --units-per-um 3.076923 \
      --out-dir ./out/pixelseq_hippo \
      --id pixelseq_hippo \
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
      --colname-transcript-x xcoord \
      --colname-transcript-y ycoord \
      --colname-transcript-feature geneName \
      --units-per-um 3.076923 \
      --out-dir /data/out/pixelseq_hippo \
      --id pixelseq_hippo \
      --width 18 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

Pixel-seq has no dedicated `run_together` profile, so the tab-separated file goes through the [generic profile](../../reference/platforms/generic.md): the `--colname-transcript-*` options name its coordinate and gene columns (each row counts as one molecule, as in the step-by-step path), and `--units-per-um` sets the pixel size. The result to deploy is `out/pixelseq_hippo/cartl/pixelseq_hippo/` (a `catalog.yaml` and the PMTiles it lists). Add `--dry-run` to print the commands without running them.

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
  include-markdown "../../../includes/includemd_vigenettes_inputformat_pixelseq.md"
%}

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="pixelseq_hippo"                # change this to reflect your dataset name
PLATFORM="pixel_seq"                    # platform information
SCALE=3.076923                          # coordinate units per micrometer

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

!!! question "How to define Scaling Factors for Pixel-Seq?"

    In [Pixel-Seq publication](https://doi.org/10.1016/j.cell.2022.10.021):
    > "Because polonies have varied sizes and shapes, to maximize the feature resolution we developed a base-calling pipeline to determine the major barcode species in each pixel (0.325 × 0.325 µm²) of gel images to construct a spatial barcode map".
    
    Accordingly, we defined scale as 1/0.325 = 3.076923

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
  include-markdown "../../../includes/includemd_vigenettes_output_pixelseq.md"
%}
