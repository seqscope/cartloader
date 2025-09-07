Generate pmtiles and web-compatible tile directories. See more details in [Reference page](../docs/reference/run_cartload2.md).

```bash
docker run -it --rm \
  -v $(pwd):/data \
  weiqiuc/cartloader:20250708b \
  run_cartload2 \
    --makefn run_cartload2.mk \
    --fic-dir /data/ficture2 \
    --out-dir /data/cartload2 \
    --id ${DATA_ID} \
    --spatula /app/cartloader/submodules/spatula/bin/spatula \
    --pmtiles /usr/local/bin/pmtiles \
    --tippecanoe /usr/local/bin/tippecanoe \
    --n-jobs ${n_jobs} \
    --threads ${n_jobs}
```