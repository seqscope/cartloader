# Data Repository Upload 

## Overview

`cartloader` offers two modules (`upload_aws` and `upload_zenodo`) to upload all generated outputs —- including rasterized SGE tiles, decoded spatial factor maps, and molecule-factor joins —- to a target data repository for sharing or deployment. It supports both AWS and Zenodo as upload backends, allowing users to choose their preferred platform.

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
    - A personal access token saved in a file, required for authentication with the Zenodo API
    - (Optional) A Zenodo deposition ID if the user prefers to upload files to an existing deposition.

    ??? "What are the Zenodo token and deposition ID, why do you need them, and how do you get them?"

        **Zenodo Token File**

        `cartloader` leverages the Zenodo API for uploading files, which provides a streamlined and efficient experience. To authenticate uploads via the API, Zenodo requires an **access token** for secure.

        To obtain a token for use with `cartloader`:

        1. Log in to [Zenodo](https://zenodo.org).
        2. Go to your [applications page](https://zenodo.org/account/settings/applications/).
        3. Click **"New Token"** and select appropriate scopes (e.g., `deposit:write`, `deposit:actions`).
        4. Copy the generated token.
        5. Save the token in a plain text file, and pass the file path to the `--token-file` option when running `cartloader`.


        **Zenodo Deposition ID**

        A **deposition ID** is a unique numeric identifier assigned to a deposition (i.e., a dataset record) you create on Zenodo. This ID specifies where your uploaded files will be stored.
        If you’ve already created a deposition, you can find its ID at the end of the URL. For example:
        ```
        https://zenodo.org/deposit/1234567
                                        ↑
                                This is the deposition ID
        ```


<!-- {%
  include-markdown "../../includes/includemd_vigenettes_upload.md"
%}
 -->
## AWS Uploads

We recommend to create create a directory within your AWS S3 bucket using a data id as directory name.

```bash
AWS_BUCKET="EXAMPLE_AWS_BUCKET"         # replace EXAMPLE_AWS_BUCKET with your actual S3 bucket name
DATA_ID="EXAMPLE_ID"                    # change EXAMPLE_ID to reflect your dataset name

cartloader upload_aws \
  --in-dir /path/to/run_cartload2/output/directory \
  --s3-dir "s3://${AWS_BUCKET}/${DATA_ID}" \
  --aws /path/to/your/aws/binary \
  --n-jobs 10
```

* `--in-dir` (str): Path to the input directory containing the cartloader compilation output
* `--s3-dir` (str): Path to the target S3 directory for uploading.
* `--aws` (str): Path to the AWS CLI binary
* `--n-jobs` (str): Number of parallel jobs
* `--catalog-yaml` (str): Path to the `catalog.yaml` file generated in `run_cartload2`. If absent, will use the `catalog.yaml` in the input directory specified by `--in-dir`.

## Zenodo Uploads
!!! warning "Define Zenodo deposition"
    You must choose one of the following modes:

      * Use an existing deposition via `--zenodo-deposition-id`
      * Create a new deposition via `--create-new-deposition`

(1) To upload files to an exist deposition ID:
```bash
zenodo_depostion_ID=DEPOSTION_ID                # Replace DEPOSTION_ID with yours

cartloader upload_zenodo \
    --in-dir /path/to/run_cartload2/output/directory \
    --upload-method catalog \
    --zenodo-token /path/to/zenodo_token.txt \
    --zenodo-deposition-id ${zenodo_depostion_ID}
```

(2) To upload files to a new deposition ID:
```bash
cartloader upload_zenodo \
  --upload-method catalog \
  --in-dir /path/to/run_cartload2/output/directory \
  --zenodo-token /path/to/zenodo_token.txt \
  --create-new-deposition \
  --title  "Title Info" \        
  --creators "Creator Name" \   
  --description "Description Info"
```

### Input Parameters
- `--in-dir` (str):Path to the input directory containing the `run_cartload2` output files.
- `--upload-method` (str, default: `all`):  Method to determine which files to upload. Options:  
    - `all`: Upload all files in `--in-dir`  
    - `catalog`: Upload files listed in a catalog YAML file  
    - `user_list`: Upload files explicitly listed via `--in-list`
- `--in-list` (list of str): Required if if using `--upload-method user_list`. Allow multiple filenames.
* `--catalog-yaml` (str): Required if `--upload-method catalog`. Path to the `catalog.yaml` file generated in `run_cartload2`. If absent, will use the `catalog.yaml` in the input directory specified by `--in-dir`.

### Zenodo Configuration
Must specify exactly one of `--zenodo-deposition-id` or `--create-new-deposition`.

* `--zenodo-token` (str): Path to your Zenodo access token file.
* `--zenodo-deposition-id` (str): A Zenodo deposition ID to upload files to. 
* `--create-new-deposition` (flag): If set, a new Zenodo deposition will be created.
* `--create-new-version` (flag): If set, a new version will be created for the provided deposition ID. This is useful to update the files for a published deposition. When enabled, `--zenodo-deposition-id` must be applied.

### Deposition Metadata
Required only if creating a new deposition.

* `--title` (str): Title for the new Zenodo deposition.
* `--upload-type` (str, default: dataset): Type of deposition. Options: dataset, software, publication, poster, presentation, image, video, lesson, other
* `--creators` (list of str): List of creators in "Lastname, Firstname" format.

### Behavior Flags:
* `--overwrite` (flag): If set, overwrite existing files in the Zenodo deposition.
* `--dry-run` (flag): If set, simulate the upload without modifying the Zenodo deposition.
<!-- * `--publish` (flag): If set, publish the deposition immediately after the upload. Recommended to review the deposition on Zenodo manually before publishing. -->
