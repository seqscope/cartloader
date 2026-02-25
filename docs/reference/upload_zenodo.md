# Upload to Zenodo

## Overview

Use `upload_zenodo` to publish `CartLoader` outputs — including PMTiles, decoded spatial factors, and the catalog — to a Zenodo deposition, either an existing one or a newly created draft. In catalog mode, files are pulled from `catalog.yaml` and grouped into cartload basics, cartload optional files (e.g., UMAP, alias), and additional basemaps (non-SGE PMTiles).

---
## Requirements

- A completed run of [`run_cartload2`](./run_cartload2.md), which produces:
    - Rasterized SGE tiles
    - (Optional) Decoded spatial factor maps
    - (Optional) Joined molecule-factor outputs
    - (Optional) Cell assets
    - (Optional) Background assets, such as histology
    - A catalog file (`catalog.yaml`) summarizing the output structure and metadata
- A Zenodo access token saved in a file (for API authentication).
- (Optional) A Zenodo deposition ID to upload to an existing deposition.

!!! info "What are the Zenodo token and deposition ID, and how do you get them?"
    **Zenodo Token File**

    `CartLoader` uses the Zenodo API for uploading files. To authenticate, Zenodo requires an access token.

    To obtain a token for use with `CartLoader`:

    1. Log in to [Zenodo](https://zenodo.org).
    2. Go to your [applications page](https://zenodo.org/account/settings/applications/).
    3. Click "New Token" and select appropriate scopes (e.g., `deposit:write`, `deposit:actions`).
    4. Copy the generated token.
    5. Save the token in a plain text file and pass the file path to `--zenodo-token` when running `CartLoader`.

    **Zenodo Deposition ID**

    A deposition ID is a unique numeric identifier assigned to a deposition (i.e., a dataset record) you create on Zenodo. This ID specifies where your uploaded files will be stored.
    If you’ve already created a deposition, you can find its ID at the end of the URL. For example:
    ```
    https://zenodo.org/deposit/1234567
                                        ↑
                                This is the deposition ID
    ```
## Actions

- Upload files to a Zenodo deposition. 
- (Optional) Create a new deposition if no `--zenodo-deposition-id` is provided or `--zenodo-deposition-id` points to a published record.
- (Optional) Fullfill or update metadata.

---
## Example Usage

=== "Create a new deposition to upload"

    ```bash
    cartloader upload_zenodo \
      --upload-method catalog \
      --in-dir /path/to/run_cartload2/output/directory \
      --zenodo-token /path/to/zenodo_token.txt \
      --title  "Title Info" \
      --creators "Creator Name" \
      --description "Description Info"
    ```

    - `--upload-method catalog`: Upload files defined in the catalog.
    - (Recommended) Provide metadata via `--title`, `--creators`, and `--description`.

=== "Upload to an existing deposition"

    ```bash
    zenodo_deposition_id=DEPOSITION_ID                # replace with your deposition ID (published IDs get a new draft automatically)

    cartloader upload_zenodo \
      --upload-method catalog \
      --in-dir /path/to/run_cartload2/output/directory \
      --zenodo-token /path/to/zenodo_token.txt \
      --zenodo-deposition-id ${zenodo_deposition_id}
    ```
    - `--upload-method catalog`: Upload files defined in the catalog.

=== "Upload explicit files"

    ```bash
    cartloader upload_zenodo \
      --upload-method files \
      --in-dir /path/to/run_cartload2/output/directory \
      --files catalog.yaml genes_all.pmtiles some_factor.pmtiles \
      --zenodo-token /path/to/zenodo_token.txt \
      --zenodo-deposition-id ${zenodo_deposition_id}
    ```
    - `--upload-method files`: Upload files defined in `--files`.

---

## Parameters

### Input Selection
- `--in-dir` (str): Path to the directory containing `run_cartload2` outputs.
- `--upload-method` (str, default: `all`): Which files to upload.

    !!! info "`--upload-method` Options"
        * `all`: Upload every file in `--in-dir`,
        * `catalog`: Upload files listed in a `catalog.yaml`,
        * `files`: Upload only the filenames provided via `--files`.

- `--files FILE [FILE ...]`: Filenames to upload (relative to `--in-dir` or absolute). Use with `--upload-method files`.
- `--catalog-yaml` (str): Path to a `catalog.yaml`. Required if `--upload-method catalog`. If omitted, uses `<in-dir>/catalog.yaml`.

### Zenodo Configuration
- `--zenodo-token` (str): Path to your Zenodo access token file.
- `--zenodo-deposition-id` (str): Zenodo deposition ID to upload files to. If the ID is for a published deposition, a new draft version is created automatically. If omitted, a new draft deposition is created. 
- `--publish` (flag): Publish the deposition automatically after upload. Recommended to publish manually after reviewing via the Zenodo web interface.

### Zenodo Metadata
- `--title` (str): Title for Zenodo deposition. Required only when creating a new deposition (i.e., if `--zenodo-deposition-id` is not provided or the existing deposition lacks these fields).
- `--upload-type` (str, default: dataset): One of: dataset, software, publication, poster, presentation, image, video, lesson, other.
- `--creators` (list of str): Creators in "Lastname, Firstname" format.

### Behavior Flags
- `--dry-run` (flag): Simulate and print operations; do not execute uploads.
- `--restart` (flag): Re-upload overlapping files instead of skipping existing ones.

---
## Output

- Files uploaded to the Zenodo deposition (existing or newly created draft).
- A `cartload.zenodo.done` marker file when `--upload-method catalog` completes the basics stage.
