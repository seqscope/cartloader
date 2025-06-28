### AWS Uploads

Copy the generated cartloader outputs to your designated AWS S3 catalog path:

```bash
cartloader upload_aws \
  --in-dir ./cartload2 \
  --s3-dir "s3://${AWS_BUCKET}/${DATA_ID}" \
  --aws ${aws} \
  --n-jobs 10
```

| Parameter       | Required  | Type   | Description                                                                 |
|-----------------|-----------|--------|-----------------------------------------------------------------------------|
| `--in-dir`      | required  | string | Path to the input directory containing the cartloader compilation output    |
| `--s3-dir`      | required  | string | Path to the target S3 directory for uploading                               |
| `--aws`         |           | string | Path to the AWS CLI binary                                                  |
| `--n-jobs`      |           | int    | Number of parallel jobs                                                     |