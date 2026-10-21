# Upload to AWS S3

!!! tip "Usually run for you"
    [`run_together`](./run_together.md#publishing) uploads each packaged sample (and, for joint runs, `multi-catalog.yaml` with the shared files) with its own `--s3-upload` stage, to `<s3-prefix>/batch=<YYYY_MM>/<collection>/<dir>`. Use `upload_aws` to upload a packaged directory directly, to any S3 path.

## Overview

Use `upload_aws` to publish `CartLoader` outputs (PMTiles, factor files, and the catalog) to Amazon S3 for sharing or deployment. It supports single‑dataset uploads and collection uploads (multiple dataset directories via `--in-list`). File lists are taken from `catalog.yaml` and split into three groups:

- **basics**: every file the catalog references, except the two groups below (includes the SGE basemaps and AI alias files recorded as `alias_ai`);
- **optional**: UMAP files (`umap.tsv`, `umap.png`, `umap.ind_png`, `umap.pmtiles`) and manual `alias` files;
- **basemaps**: non-SGE basemaps, such as histology images.

!!! question
    - [Why did my AWS upload fail, and how can I fix it?](../faq/aws_upload_failures.md)

---
## Requirements

- A packaged directory with a `catalog.yaml`, from [`run_cartload2`](./run_cartload2.md) (or each sample directory of [`run_cartload2_multi`](./run_cartload2_multi.md)). Every file the catalog references must exist locally.
- AWS CLI installed and configured (e.g., `aws configure`).

---
## Actions

Writes a Makefile (`--makefn`, inside `--in-dir`) with one upload target per group and runs it:

- basics: all basic files plus `catalog.yaml` in one target, flagged by `<in-dir>/cartload.aws.done`;
- optional and basemap files: one target per file, each followed by a re-upload of `catalog.yaml`, flagged by `<file>.aws.done`.

`catalog.yaml` is always uploaded as `<s3-dir>/catalog.yaml`. Existing flags make reruns skip completed uploads; use `--restart` to upload again.

---
## Example Usage

=== "Single Dataset Upload"
    ```bash
    AWS_DIR="s3://example/aws/prefix"   # S3 prefix that will contain the dataset
    DATA_ID="EXAMPLE_ID"                # dataset identifier

    cartloader upload_aws \
      --in-dir /path/to/cartload2/output \
      --s3-dir "${AWS_DIR}/${DATA_ID}" \
      --n-jobs 10
    ```

=== "Collection Upload (Multiple Datasets)"

    ```bash
    AWS_DIR="s3://example/aws/prefix"        # S3 prefix that will contain the collection
    COLLECTION_ID="EXAMPLE_COLLECTION_ID"    # collection identifier
    IN_LIST=/path/to/dirs.tsv                # one sub-directory name per line; no header
    cartloader upload_aws \
      --in-dir /path/to/parent_cartload_outputs \
      --s3-dir "${AWS_DIR}/${COLLECTION_ID}" \
      --in-list ${IN_LIST} \
      --n-jobs 10
    ```

    !!! info "Collection Structure"

        **Input Collection**

          * `--in-dir` is a parent directory containing one subdirectory per dataset; `--in-list` names those subdirectories (first column).
          * Each subdirectory must include its own `catalog.yaml`.
          * For a `run_cartload2_multi` output, the subdirectories are named `<multi_id>-<sample_id>`; list those names. Collection mode uploads only the subdirectories, not the root `multi-catalog.yaml` or the shared factor files next to it.

        ```
        /path/to/parent_cartload_outputs/   # defined by `--in-dir`
        ├── SAMPLE_001      # listed in `--in-list`
        │   ├── catalog.yaml
        │   ├── ...         # SGE assets (e.g., genes_all.pmtiles)
        │   ├── ...         # FICTURE assets (e.g., *-model.tsv, *-pseudobulk.tsv.gz, *-bulk-de.tsv, *-info.tsv, *.pmtiles, *-rgb.tsv)
        │   ├── ...         # Background assets (e.g., from `cartloader import_image`)
        │   └── ...         # Cell assets (e.g., from `import_*_cell`)
        ├── SAMPLE_002
        │   └── ...
        └── SAMPLE_003
           └── ...
        ```

        **Output S3 Collection**

        * Each dataset uploads to `<s3-dir>/<s3_id>/`, where `<s3_id>` is the catalog's `id`, or, without one, the listed name lowercased with underscores replaced by hyphens.

        ```
        s3://example/aws/prefix/EXAMPLE_COLLECTION_ID/   # defined by `--s3-dir`
        ├── sample-001      # catalog `id`, else lowercased name with underscores → hyphens
        │   ├── catalog.yaml
        │   └── ...
        ├── sample-002
        │   └── ...
        └── sample-003
           └── ...
        ```

---

## Parameters

### Input/Output
- `--in-dir` (str, required): Input directory (single) or parent directory with per-dataset subdirectories (collection).
- `--s3-dir` (str, required): S3 destination (single) or parent prefix for per-dataset outputs (collection).
- `--catalog-yaml` (str): Catalog to read (single mode only; default: `<in-dir>/catalog.yaml`).

### Upload Selection
- `--upload-basics-only` (flag): Upload only the basic files and the catalog.
- `--upload-optional-only` (flag): Upload only the optional files (UMAP, alias).
- `--upload-basemap-only` (flag): Upload only the additional basemaps.

Use at most one of these flags: combining two or more uploads nothing.

### Collection Parameters
- `--in-list` (str): TSV of subdirectory names under `--in-dir` (first column; no header); enables collection mode. Cannot be combined with `--catalog-yaml`.

### AWS Configuration
- `--aws` (str): AWS CLI executable (default: `aws`).
- `--profile` (str): AWS CLI profile passed to every `aws s3 cp` (default: the CLI's default profile).

### Run Parameters
- `--dry-run` (flag): Write the Makefile and print the commands (`make -n`) without uploading.
- `--restart` (flag): Ignore existing `.aws.done` flags and upload everything again.
- `--n-jobs` (int): Number of parallel jobs (default: `2`).
- `--makefn` (str): Makefile name, written inside `--in-dir` (default: `upload_aws.mk`).

---
## Output

- Files uploaded under `--s3-dir` (per-dataset subdirectories in collection mode).
- `<in-dir>/upload_aws.mk` (always written, also in dry-run mode).
- Flag files: `cartload.aws.done` in each dataset directory and `<file>.aws.done` next to each uploaded optional or basemap file.
