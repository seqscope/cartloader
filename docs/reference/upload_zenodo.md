# Upload to Zenodo

## Overview

Use `upload_zenodo` to publish CartLoader outputs — including PMTiles, decoded spatial factors, and the catalog — to a Zenodo deposition, either an existing one or a newly created draft.

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

!!! note "What are the Zenodo token and deposition ID, and how do you get them?"
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

## Example Usage

1) Upload files to an existing deposition ID:

```bash
zenodo_deposition_id=DEPOSITION_ID                # replace with your deposition ID

cartloader upload_zenodo \
  --in-dir /path/to/run_cartload2/output/directory \
  --upload-method catalog \
  --zenodo-token /path/to/zenodo_token.txt \
  --zenodo-deposition-id ${zenodo_deposition_id}
```

2) Create a new deposition and upload files (supply metadata):

```bash
cartloader upload_zenodo \
  --upload-method catalog \
  --in-dir /path/to/run_cartload2/output/directory \
  --zenodo-token /path/to/zenodo_token.txt \
  --title  "Title Info" \        
  --creators "Creator Name" \   
  --description "Description Info"
```

## Parameters

### Input Selection
- `--in-dir` (str): Path to the directory containing `run_cartload2` outputs.
- `--upload-method` (str, default: `all`): Choose one of the following input modes to define files to upload:
    - `all`: Upload every file in `--in-dir`.
    - `catalog`: Upload files listed in a catalog YAML file.
    - `user_list`: Upload only files specified via `--in-list`.
- `--in-list` (list of str): List of files. Required if using `--upload-method user_list`.
- `--catalog-yaml` (str): Path to a catalog.yaml file. Required if `--upload-method catalog`. If omitted, uses `<in-dir>/catalog.yaml`.

### Zenodo Configuration
- `--zenodo-token` (str): Path to your Zenodo access token file.
- `--zenodo-deposition-id` (str): Zenodo deposition ID to upload files to. If the ID is for a published deposition, a new draft version is created automatically. If omitted, a new draft deposition is created. 
- `--title` (str): Title for the new Zenodo deposition. Required only when creating a new deposition (i.e., if `--zenodo-deposition-id` is not provided or the existing deposition lacks these fields).
- `--upload-type` (str, default: dataset): One of: dataset, software, publication, poster, presentation, image, video, lesson, other.
- `--creators` (list of str): Creators in "Lastname, Firstname" format.
- `--publish` (flag): Publish the deposition automatically after upload. Recommended to publish manually after reviewing via the Zenodo web interface.

### Behavior Flags
- `--dry-run` (flag): Simulate the upload without transferring files.
- `--restart` (flag): Re-upload files regardless of existing ones.

