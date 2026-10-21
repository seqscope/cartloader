=== "AWS Uploads"

    Upload the generated `CartLoader` outputs to your designated AWS S3 directory. The container needs your AWS credentials, so mount `~/.aws` read-only:
    ```bash
    # AWS S3 target location
    S3_DIR=s3://your-bucket/path/to/dir     # recommended: end with DATA_ID, e.g. s3://your-bucket/test-data/${DATA_ID}

    docker run --rm \
      -v $(pwd)/out:/out \
      -v ${HOME}/.aws:/root/.aws:ro \
      hyunminkang/cartloader:latest \
      upload_aws \
        --in-dir /out/cartload2 \
        --s3-dir "${S3_DIR}" \
        --n-jobs ${n_jobs}
    ```

    | Parameter   | Required | Type   | Description |
    |-------------|----------|--------|-------------|
    | `--in-dir`  | required | string | Directory with the `CartLoader` packaging output (containing `catalog.yaml`); a container path |
    | `--s3-dir`  | required | string | Target S3 directory |
    | `--profile` |          | string | AWS CLI profile to use (from the mounted `~/.aws`) |
    | `--n-jobs`  |          | int    | Number of parallel jobs (default: `2`) |


=== "Zenodo Uploads"

    Upload the generated `CartLoader` outputs to a new Zenodo deposition or an existing one. Mount the token file into the container and pass its container path:

    ```bash
    zenodo_token=/path/to/zenodo/token/file    # file holding your Zenodo access token (on the host)

    docker run --rm \
      -v $(pwd)/out:/out \
      -v ${zenodo_token}:/root/zenodo_token:ro \
      hyunminkang/cartloader:latest \
      upload_zenodo \
        --in-dir /out/cartload2 \
        --upload-method catalog \
        --zenodo-token /root/zenodo_token \
        --title "Your Title" \
        --creators "Lastname, Firstname" \
        --description "This is an example description"
    ```

    | Parameter                | Required              | Type        | Description |
    |--------------------------|-----------------------|-------------|-------------|
    | `--in-dir`               | required              | string      | Directory with the `CartLoader` packaging output; a container path |
    | `--zenodo-token`         | required              | string      | Container path of the mounted token file |
    | `--upload-method`        |                       | string      | Which files to upload: `all` (every file in `--in-dir`; default), `catalog` (the files listed in `catalog.yaml`), or `files` (the names given with `--files`) |
    | `--zenodo-deposition-id` |                       | string      | Existing deposition to upload into (default: create a new one) |
    | `--title`                | required <sup>1</sup> | string      | Title of the deposition |
    | `--creators`             | required <sup>1</sup> | list of str | Creators in "Lastname, Firstname" format, each quoted |

    <sub><sup>1</sup>: Required when creating a new deposition, or when the existing one lacks it.</sub>
