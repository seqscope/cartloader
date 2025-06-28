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
    - A Zenodo deposition ID

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


## Example Usage

<!-- {%
  include-markdown "../../includes/includemd_vigenettes_upload2aws.md"
%}
 -->
### AWS Uploads

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

### Zenodo Uploads

```bash
zenodo_depostion_ID=DEPOSTION_ID                # Replace DEPOSTION_ID with yours

cartloader upload_zenodo \
    --zenodo-id ${zenodo_depostion_ID}
    --zenodo-token /path/to/zenodo_token.txt \
    --in-dir /path/to/run_cartload2/output/directory \
    --upload-method catalog
```

* `--in-dir` (str): Path to the input directory containing the cartloader compilation output
* `--zenodo-id` (str): Zenodo Deposition ID
* `--zenodo-token` (str): Path to the Zenodo token file
* `--upload-method` (str): Specifiy `catalog` to use the catalog file to locate files to be uploaded
* `--catalog-yaml` (str): Path to the `catalog.yaml` file generated in `run_cartload2`. If absent, will use the `catalog.yaml` in the input directory specified by `--in-dir`.
