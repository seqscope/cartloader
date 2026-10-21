```bash
# ====
# Replace the placeholder with the actual path on your system.
# ====
work_dir=/path/to/work/directory        # directory holding the input data; mounted at /data in the container
cd $work_dir

n_jobs=10                               # number of parallel jobs (also used as threads per job below)

docker pull hyunminkang/cartloader:latest
```

!!! info "Using the Docker image"
    - The image's entry point is `cartloader`: whatever follows the image name is a `cartloader` subcommand and its options (`docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest <command> [options]`).
    - Paths passed to the command are **container** paths: with `-v ${work_dir}:/data`, a host file `${work_dir}/raw/x.csv` is `/data/raw/x.csv`.
    - Every tool is built into the image, so no tool-path options are needed.
    - Uploads need credentials inside the container: mount them read-only, e.g. `-v ${HOME}/.aws:/root/.aws:ro` for AWS.
    - `latest` is the most recent release; for reproducible runs use a dated tag (e.g. `hyunminkang/cartloader:20261007-1432`).
