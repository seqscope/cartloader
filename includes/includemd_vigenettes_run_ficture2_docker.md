Compute spatial factors using `punkst` (FICTURE2). See more details in [Reference page](../docs/reference/run_ficture2.md).

```bash
docker run -it --rm \
  -v $(pwd):/data \
  weiqiuc/cartloader:${docker_tag} \
  run_ficture2 \
    --makefn run_ficture2.mk \
    --main \
    --in-transcript /app/data/sge/transcripts.unsorted.tsv.gz \
    --in-feature /app/data/sge/feature.clean.tsv.gz \
    --in-minmax /app/data/sge/coordinate_minmax.tsv \
    --cmap-file ${cmap} \
    --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \
    --out-dir /data/ficture2 \
    --width ${train_width} \
    --n-factor ${n_factor} \
    --ficture2 ${punkst} \
    --spatula ${spatula} \
    --n-jobs ${n_jobs} \
    --threads ${n_jobs}
```
