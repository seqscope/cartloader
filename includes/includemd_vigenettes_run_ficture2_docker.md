Compute spatial factors using `punkst` (FICTURE2 mode). See more details in [Reference page](../docs/reference/run_ficture2.md).

```bash
docker run -it --rm \
  -v $(pwd):/data \
  weiqiuc/cartloader:20250708b \
  run_ficture2 \
    --makefn run_ficture2.mk \
    --main \
    --in-transcript /app/data/sge/transcripts.unsorted.tsv.gz \
    --in-feature /app/data/sge/feature.clean.tsv.gz \
    --in-minmax /app/data/sge/coordinate_minmax.tsv \
    --cmap-file /app/cartloader/assets/fixed_color_map_256.tsv \
    --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \
    --out-dir /data/ficture2 \
    --width ${train_width} \
    --n-factor ${n_factor} \
    --ficture2 /app/cartloader/submodules/punkst \
    --spatula /app/cartloader/submodules/spatula/bin/spatula \
    --n-jobs ${n_jobs} \
    --threads ${n_jobs}
```
