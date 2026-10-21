Compute spatial factors using `punkst` (FICTURE2). See more details on the [Reference page](../docs/reference/run_ficture2.md).

This example reads the Seq-Scope starter SGE bundled in the image (`/app/data/sge/`) and writes to `./out` on the host, mounted at `/out`. To use your own SGE, mount its directory too and change the `--in-*` paths to the container paths.

```bash
mkdir -p out
docker run --rm \
  -v $(pwd)/out:/out \
  hyunminkang/cartloader:latest \
  run_ficture2 \
    --makefn run_ficture2.mk \
    --main \
    --in-transcript /app/data/sge/transcripts.unsorted.tsv.gz \
    --in-feature /app/data/sge/feature.clean.tsv.gz \
    --in-minmax /app/data/sge/coordinate_minmax.tsv \
    --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \
    --out-dir /out/ficture2 \
    --width ${train_width} \
    --n-factor ${n_factor} \
    --n-jobs ${n_jobs} \
    --threads ${n_jobs}
```
