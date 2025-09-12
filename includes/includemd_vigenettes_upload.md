<!--section1-start-->

!!! info "Choose a data repository to host/share your output"

    `CartLoader` supports two upload options (`AWS` and `Zenodo`) for storing PMTiles of SGE and spatial factors in a data repository.

    **Choose the one that best suits your needs.**

### AWS Uploads

Upload the generated CartLoader outputs to your designated AWS S3 directory:
<!--section1-end-->

```bash
# AWS S3 target location for cartostore
S3_DIR=/s3/path/to/s3/dir              # Recommend to use DATA_ID as directory name, such as s3://bucket_name/xenium-v1-humanlung-cancer-ffpe

cartloader upload_aws \
  --in-dir ./cartload2 \
  --s3-dir "${S3_DIR}" \
  --aws ${aws} \
  --n-jobs ${n_jobs}
```

| Parameter       | Required  | Type   | Description                                                                 |
|-----------------|-----------|--------|-----------------------------------------------------------------------------|
| `--in-dir`      | required  | string | Path to the input directory containing the CartLoader asset packaging output    |
| `--s3-dir`      | required  | string | Path to the target S3 directory for uploading                               |
| `--aws`         |           | string | Path to the AWS CLI binary                                                  |
| `--n-jobs`      |           | int    | Number of parallel jobs                                                     |

<!--section2-start-->
### Zenodo Uploads

Upload the generated CartLoader outputs to your designated Zenodo deposition or a new deposition.
<!--section2-end-->

```bash
zenodo_token=/path/to/zenodo/token/file    # replace /path/to/zenodo/token/file by path to your zenodo token file

cartloader upload_zenodo \
  --in-dir ./cartload2 \
  --upload-method catalog \
  --zenodo-token $zenodo_token \
  --title  "Your Title" \
  --creators "Your Name" \
  --description "This is an example description"
```
<!--section3-start-->
| Parameter         | Required | Type        | Description                                                                                                                                                                                                         |
|-------------------|----------|-------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--in-dir`        | required | string      | Path to the input directory containing the CartLoader asset packaging output                                                                                                                                            |
| `--upload-method` | required | string      | Method to determine which files to upload. Options: `all` to upload all files in `--in-dir`; `catalog` to upload files listed in a catalog YAML file; `user_list` to upload files explicitly listed via `--in-list` |
| `--catalog-yaml`  |          | string      | Required if `--upload-method catalog`. Path to `catalog.yaml` generated in `run_cartload2`. If absent, uses the catalog in the input directory specified by `--in-dir`.                                             |
| `--zenodo-token` | required | string      | Path to your Zenodo access token file                                                                                                                                                                               |
| `--title`         | required | string      | Required when creating a new deposition (i.e., if `--zenodo-deposition-id` is omitted). Title for the new Zenodo deposition.                                                                                        |
| `--creators`      | required | list of str | List of creators in "Lastname, Firstname" format.                                                                                                                                                                   |
<!--section3-end-->
