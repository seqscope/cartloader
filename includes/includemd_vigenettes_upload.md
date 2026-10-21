<!--section1-start-->
!!! info "Choose a data repository to host/share your output"

    `CartLoader` supports two upload options (`AWS` and `Zenodo`) for storing PMTiles of SGE and spatial factors in a data repository.

    **Choose the one that best suits your needs.**
<!--section1-end-->

=== "AWS Uploads"
    Upload the generated `CartLoader` outputs to your designated AWS S3 directory. The `aws` CLI must be configured with credentials that can write there (see [AWS Upload Failures](../docs/faq/aws_upload_failures.md)).

    ```bash
    # AWS S3 target location
    S3_DIR=s3://your-bucket/path/to/dir     # recommended: end with DATA_ID, e.g. s3://your-bucket/test-data/${DATA_ID}

    cartloader upload_aws \
      --in-dir ./cartload2 \
      --s3-dir "${S3_DIR}" \
      --n-jobs ${n_jobs}
    ```

    | Parameter   | Required | Type   | Description |
    |-------------|----------|--------|-------------|
    | `--in-dir`  | required | string | Directory with the `CartLoader` packaging output (containing `catalog.yaml`) |
    | `--s3-dir`  | required | string | Target S3 directory |
    | `--profile` |          | string | AWS CLI profile to use |
    | `--aws`     |          | string | Path to the AWS CLI (default: `aws`) |
    | `--n-jobs`  |          | int    | Number of parallel jobs (default: `2`) |

=== "Zenodo Uploads"
    Upload the generated `CartLoader` outputs to a new Zenodo deposition or an existing one (see [Zenodo Setup](../docs/faq/zenodo_setup.md) for the access token).

    ```bash
    zenodo_token=/path/to/zenodo/token/file    # file holding your Zenodo access token

    cartloader upload_zenodo \
      --in-dir ./cartload2 \
      --upload-method catalog \
      --zenodo-token ${zenodo_token} \
      --title "Your Title" \
      --creators "Lastname, Firstname" \
      --description "This is an example description"
    ```

    | Parameter                | Required              | Type        | Description |
    |--------------------------|-----------------------|-------------|-------------|
    | `--in-dir`               | required              | string      | Directory with the `CartLoader` packaging output |
    | `--zenodo-token`         | required              | string      | Path to the file holding your Zenodo access token |
    | `--upload-method`        |                       | string      | Which files to upload: `all` (every file in `--in-dir`; default), `catalog` (the files listed in `catalog.yaml`), or `files` (the names given with `--files`) |
    | `--catalog-yaml`         |                       | string      | Catalog to read with `--upload-method catalog` (default: `<in-dir>/catalog.yaml`) |
    | `--zenodo-deposition-id` |                       | string      | Existing deposition to upload into (default: create a new one) |
    | `--title`                | required <sup>1</sup> | string      | Title of the deposition |
    | `--creators`             | required <sup>1</sup> | list of str | Creators in "Lastname, Firstname" format, each quoted |
    | `--description`          |                       | string      | Description of the deposition |

    <sub><sup>1</sup>: Required when creating a new deposition, or when the existing one lacks it.</sub>
