```bash
# AWS S3 target location for cartostore
AWS_BUCKET="EXAMPLE_AWS_BUCKET"         # replace EXAMPLE_AWS_BUCKET with your actual S3 bucket name

docker run -it --rm \
  -v $(pwd):/data \
  weiqiuc/cartloader:${docker_tag} \
  upload_aws \
    --in-dir /data/cartload2 \
    --s3-dir "s3://${AWS_BUCKET}/${DATA_ID}" \
    --aws ${aws} \
    --n-jobs ${n_jobs}

```