```bash
docker run -it --rm \
  -v $(pwd):/data \
  weiqiuc/cartloader:20250708b \
  upload_aws \
    --in-dir /data/cartload2 \
    --s3-dir "s3://${AWS_BUCKET}/${DATA_ID}" \
    --aws /usr/local/bin/aws \
    --n-jobs ${n_jobs}

```