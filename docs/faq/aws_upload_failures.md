# Why did my AWS upload fail, and how can I fix it?

If `cartloader upload_aws` fails, it is usually one of the issues below.

## 1) AWS credentials are missing or invalid

Common messages:

- `Unable to locate credentials`
- `The security token included in the request is invalid`
- `ExpiredToken`

Checks:

```bash
aws sts get-caller-identity
aws configure list
```

Fix:

- Configure credentials (`aws configure`) or set `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, and (if needed) `AWS_SESSION_TOKEN`.
- If you use temporary credentials, refresh them before running upload.

## 2) Permission denied on bucket/prefix

Common messages:

- `AccessDenied`
- `An error occurred (AccessDenied) when calling ...`

Checks:

```bash
aws s3 ls s3://your-bucket/
aws s3 ls s3://your-bucket/your-prefix/
```

Fix:

- Ensure IAM policy allows at least `s3:ListBucket`, `s3:PutObject`, and `s3:GetObject` on the target bucket/prefix.
- Confirm bucket policy does not block your principal.

## 3) Wrong S3 destination path

Common issue:

- Upload appears to "succeed" but files are not where you expected.

Fix:

- Verify `--s3-dir` is a full S3 URI like `s3://bucket/path`.
- In collection mode (`--in-list`), remember each sample is uploaded to a subdirectory under `--s3-dir`.

## 4) Input directory/catalog problems

Common messages:

- Missing `catalog.yaml`
- Missing expected PMTiles/assets

Fix:

- Confirm `--in-dir` points to a completed `run_cartload2` (or collection parent for `run_cartload2_multi`).
- In single mode, verify `catalog.yaml` exists (or pass `--catalog-yaml` explicitly).

## 5) Docker run cannot see AWS credentials

If local `aws` works but Docker upload fails, credentials may not be available inside the container.

Fix:

- Mount/pass AWS credentials into the container (for example, `~/.aws`) or pass required `AWS_*` env vars.
- Confirm the `aws` binary inside container is available (or set `--aws`).

## 6) Region, endpoint, or network issues

Common messages:

- `Could not connect to the endpoint URL`
- timeout/name resolution errors

Fix:

- Verify internet/VPN/proxy settings.
- Check bucket region and AWS CLI region configuration.

## Minimal preflight before upload

```bash
aws sts get-caller-identity
aws s3 ls s3://your-bucket/
```

If both pass, retry:

```bash
cartloader upload_aws --in-dir /path/to/cartload2 --s3-dir s3://your-bucket/your-prefix
```

See also:

- [AWS upload reference](../reference/upload_aws.md)
