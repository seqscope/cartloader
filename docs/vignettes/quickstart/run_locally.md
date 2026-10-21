# Quick Start: Run Locally

This quick start runs [`run_together`](../../reference/run_together.md) with a local installation of `CartLoader`, on the same two examples as the [Docker quick start](./run_in_docker.md).

!!! info "Use Cases"
    Use a local installation to follow the latest code on GitHub, to modify `CartLoader`, or when Docker is not available.

!!! warning "Requirements"
    - `CartLoader` and its tools installed and built as described in the [Installation Guide](../../installation.md), including `bash build.sh` in `submodules/`.
    - For the second example: about 10 GB of free disk (a 1.5 GB download plus about 8 GB of outputs); it took about 11 minutes on a 12-core laptop with `--threads 4 -j 4`.

!!! question
    - [I am not sure whether to run with Docker or run locally. Which one should I choose?](../../faq/choose_run_mode.md)

---

## Check the Installation

Activate the environment you installed `CartLoader` into, then check that every dependency is found:

```bash
python installation/check_dependencies.py    # run from the cartloader repository
```

`run_together` uses the tools built under `submodules/` (`punkst`, `spatula`, `pmpoint`, `tippecanoe`, `pmtiles`, `geotiff2pmtiles`), so you never pass tool paths to it.

!!! warning "After updating `CartLoader`, reinstall and rebuild"
    The pipeline expects the Python packages and tool versions pinned by the repository. After `git pull`, reinstall the Python dependencies and rebuild the submodules; otherwise a step may fail with a missing-package error (for example `requires the 'imagecodecs' package`), an "unknown option" error, or a usage message:
    ```bash
    uv pip install -e ".[ai]"                 # or: python -m pip install -e ".[ai]"
    git submodule update --init --recursive
    cd submodules && bash build.sh
    ```

---

## Example 1: Mouse Hippocampus (Seq-Scope)

### Input

Download the starter dataset, a mouse hippocampus SGE already converted to the `FICTURE`-compatible transcript format by [`sge_convert`](../../reference/sge_convert.md):

```bash
work_dir=/path/to/work/directory
mkdir -p ${work_dir} && cd ${work_dir}

wget https://zenodo.org/records/17953582/files/seqscope_starter.std.tar.gz
tar -xzf seqscope_starter.std.tar.gz      # creates ./sge/
```

{%
  include-markdown "../../../includes/includemd_vigenettes_sgeformat.md"
%}

### Run

```bash
cartloader run_together \
    --platform seqscope \
    --in-transcript ./sge/transcripts.unsorted.tsv.gz \
    --out-dir ./out/seqscope_hippo \
    --id seqscope_hippo \
    --width 18 \
    --n-factor 6,12 \
    --threads 4 -j 4
```

| Option | Meaning |
|---|---|
| `--platform seqscope` | Use the Seq-Scope [platform profile](../../reference/supported_platforms.md) for defaults. |
| `--in-transcript` | An already-converted transcript TSV, so the SGE conversion step is skipped. |
| `--out-dir`, `--id` | Output directory and the sample ID used to name the packaged assets. |
| `--width`, `--n-factor` | `FICTURE` hexagon width (µm) and the number(s) of factors; one model per combination. |
| `--threads`, `-j` | Threads per job and number of parallel jobs. |

The deployable result is `out/seqscope_hippo/cartl/seqscope_hippo/`: a `catalog.yaml` plus the PMTiles it references. See [Outputs](./run_in_docker.md#outputs) for the full layout.

---

## Example 2: Xenium Human Kidney (full pipeline)

This example starts from raw Xenium output, so every stage runs: SGE conversion, `FICTURE`, two cell analyses (`CartLoader`'s own, from the Xenium cell boundaries, and an import of the Xenium Ranger clusters), the DAPI image, and packaging.

### Input

A human kidney section from [GSE264334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264334), repackaged in the standard Xenium output layout and hosted on the [AWS Open Data Registry](https://registry.opendata.aws/) (no AWS account needed):

```bash
# with the AWS CLI
aws s3 cp --no-sign-request \
  "s3://cartostore/data/batch=2026_05/xenium-geo-public-dataset-collection/xenium-human-kidney-igan-bull2024-20260516/gse264334-reformatted.tar.gz" .
# or with curl
curl -O "https://cartostore.s3.amazonaws.com/data/batch=2026_05/xenium-geo-public-dataset-collection/xenium-human-kidney-igan-bull2024-20260516/gse264334-reformatted.tar.gz"

tar -xzf gse264334-reformatted.tar.gz    # creates gse264334-reformatted/
```

### Run

```bash
cartloader run_together \
    --platform 10x_xenium \
    --in-dir ./gse264334-reformatted \
    --out-dir ./out/gse264334 \
    --id gse264334 \
    --bin-count 50 \
    --threads 4 -j 4
```

Without `--width`/`--n-factor`, the Xenium profile's defaults apply (12 µm hexagons; 12, 24, and 48 factors). `--bin-count 50` groups genes into 50 bins for the point layers instead of the default 500, which makes packaging faster for a quick test. The deployable result is `out/gse264334/cartl/gse264334/`.

!!! tip "Interrupted or failed run?"
    Rerun the same command: finished steps are skipped. Add `--restart` to rebuild everything. See [Resume or Restart](../../faq/restart_or_resume.md).

---

## Next Steps

- Run your own data: [platform tutorials](../intro.md#platform-tutorials) and [Specifying Inputs](../../reference/run_together_inputs.md).
- Analyze several samples jointly: [Multiple Samples with run_together](../pipelines/run_together_multi.md).
- Publish to S3 and annotate factors with an LLM: [Publishing](../../reference/run_together.md#publishing).
- Run the stages one at a time: [Step-by-step](../step_by_step.md).
