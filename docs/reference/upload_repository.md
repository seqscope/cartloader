# Data Repository Upload

## Overview

`CartLoader` offers two modules (`upload_aws` and `upload_zenodo`) to upload all generated outputs — including rasterized SGE tiles, decoded spatial factor maps, and molecule–factor joins — to a target data repository for sharing or deployment. It supports both AWS and Zenodo as upload backends, allowing you to choose your preferred platform.

<!-- The uploaded package includes all FICTURE analysis outputs, PMTiles for CartoScope visualization, and catalog metadata for downstream access and integration. -->

## Requirements

- A completed run of [`run_cartload2`](./run_cartload2.md), which produces:
    - Rasterized SGE tiles
    - Decoded spatial factor maps
    - Joined molecule-factor outputs
    - A catalog file (`catalog.yaml`) summarizing the output structure and metadata
- **For AWS uploads**:
    - AWS CLI installed and configured (e.g., via `aws configure`)
- **For Zenodo uploads**:
    - A personal access token saved in a file (required to authenticate with the Zenodo API)
    - (Optional) A Zenodo deposition ID if you prefer to upload files to an existing deposition

    ??? "What are the Zenodo token and deposition ID, and how do you get them?"

        **Zenodo Token File**

        `CartLoader` uses the Zenodo API for uploading files. To authenticate, Zenodo requires an access token.

        To obtain a token for use with `CartLoader`:

        1. Log in to [Zenodo](https://zenodo.org).
        2. Go to your [applications page](https://zenodo.org/account/settings/applications/).
        3. Click **"New Token"** and select appropriate scopes (e.g., `deposit:write`, `deposit:actions`).
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

{%
  include-markdown "../../includes/includemd_vigenettes_upload.md" start="<!--section1-start-->" end="<!--section1-end-->" preserve-includer-indent=false
%}

```bash
AWS_BUCKET="EXAMPLE_AWS_BUCKET"         # replace EXAMPLE_AWS_BUCKET with your actual S3 bucket name
DATA_ID="EXAMPLE_ID"                    # change EXAMPLE_ID to reflect your dataset name

cartloader upload_aws \
  --in-dir /path/to/run_cartload2/output/directory \
  --s3-dir "s3://${AWS_BUCKET}/${DATA_ID}" \
  --aws /path/to/your/aws/binary \
  --n-jobs 10
```

{%
  include-markdown "../../includes/includemd_vigenettes_upload.md" start="<!--section2-start-->" end="<!--section2-end-->" preserve-includer-indent=false
%}


(1) Upload files to an existing deposition ID:
```bash
zenodo_deposition_id=DEPOSITION_ID                # Replace DEPOSITION_ID with your own

cartloader upload_zenodo \
    --in-dir /path/to/run_cartload2/output/directory \
    --upload-method catalog \
    --zenodo-token /path/to/zenodo_token.txt \
    --zenodo-deposition-id ${zenodo_deposition_id}
```

(2) Create a new deposition and upload files (supply metadata):
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

### Input Parameters
- `--in-dir` (str): Path to the input directory containing the `run_cartload2` output files.
- `--upload-method` (str, default: `all`): Which files to upload:
    - `all`: Upload every file in `--in-dir`.
    - `catalog`: Upload files listed in a catalog YAML file.
    - `user_list`: Upload only files specified via `--in-list`.
- `--in-list` (list of str): Required if using `--upload-method user_list`. Allows multiple filenames.
- `--catalog-yaml` (str): Required if `--upload-method catalog`. Path to `catalog.yaml` from `run_cartload2`. If absent, uses `<in-dir>/catalog.yaml`.

### Zenodo Configuration
* `--s3-dir` (str): Path to your AWS directory.

### Zenodo Configuration
* `--zenodo-token` (str): Path to your Zenodo access token file.
* `--zenodo-deposition-id` (str): Zenodo deposition ID to upload files to. If the ID is for a published deposition, a new draft version is created automatically. If omitted, a new draft deposition is created. 
* `--title` (str): Title for the new Zenodo deposition. Required only when creating a new deposition (i.e., if `--zenodo-deposition-id` is not provided or the existing deposition lacks these fields).
* `--upload-type` (str, default: dataset): Type of deposition. Options: dataset, software, publication, poster, presentation, image, video, lesson, other
* `--creators` (list of str): List of creators in "Lastname, Firstname" format.
* `--publish` (flag): Publish the deposition automatically after upload. Recommended to publish manually after reviewing via the Zenodo web interface.

### Behavior Flags
* `--dry-run` (flag): Simulate the upload without modifying the Zenodo deposition.
* `--restart` (flag): Restart all uploading process regardless of existing files.

{%
  include-markdown "../../includes/includemd_vigenettes_upload.md" start="<!--section3-start-->" end="<!--section3-end-->" preserve-includer-indent=false
%}
