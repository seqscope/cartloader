# Upload to AWS S3

## Overview

Use `upload_aws` to publish CartLoader outputs — including PMTiles, decoded spatial factors, and the catalog — to an S3 location for sharing or deployment.

## Requirements

- A completed run of [`run_cartload2`](./run_cartload2.md), which produces:
    - Rasterized SGE tiles
    - (Optional) Decoded spatial factor maps
    - (Optional) Joined molecule-factor outputs
    - (Optional) Cell assets
    - (Optional) Background assets, such as histology
    - A catalog file (`catalog.yaml`) summarizing the output structure and metadata
- AWS CLI installed and configured (e.g., `aws configure`).

## Example Usage

```bash
AWS_BUCKET="EXAMPLE_AWS_BUCKET"         # replace with your S3 bucket name
DATA_ID="EXAMPLE_ID"                    # dataset identifier

cartloader upload_aws \
  --in-dir /path/to/run_cartload2/output/directory \
  --s3-dir "s3://${AWS_BUCKET}/${DATA_ID}" \
  --aws /path/to/aws \
  --n-jobs 10
```

## Parameters

### Input/Output
- `--in-dir` (str): Path to the directory containing `run_cartload2` outputs.

### AWS Configuration
- `--s3-dir` (str): Target S3 directory (e.g., `s3://bucket/prefix`).
- `--aws` (str): Path to AWS CLI binary (default: `aws`).
- `--n-jobs` (int): Number of parallel upload jobs (default: 1).

### Behavior Flags
- `--dry-run` (flag): Simulate the upload without writing to S3.
- `--restart` (flag): Ignore existing uploaded files and re-upload.

