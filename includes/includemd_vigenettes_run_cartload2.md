Generate pmtiles and web-compatible tile directories. See more details in [Reference page](../docs/reference/run_cartload2.md).

<!-- ```bash
cartloader run_cartload2 \
    --makefn run_cartload2.mk \         # (optional) file name of the output make file
    --fic-dir ./ficture2 \              # path to input directory containing FICTURE2 output
    --out-dir ./cartload2 \             # path to output directory for PMTiles and web tiles
    --id ${DATA_ID} \                   # dataset ID used for naming outputs and metadata
    --spatula ${spatula} \              # (optional) path to the spatula binary
    --pmtiles ${pmtiles} \              # (optional) path to the pmtiles binary
    --tippecanoe ${tippecanoe} \        # (optional) path to the tippecanoe binary
    --n-jobs 10 \                       # (optional) number of parallel jobs
    --threads 10                        # (optional) number of threads per job
``` -->

```bash
cartloader run_cartload2 \
  --makefn run_cartload2.mk \
  --fic-dir ./ficture2 \
  --out-dir ./cartload2 \
  --id ${DATA_ID} \
  --spatula ${spatula} \
  --pmtiles ${pmtiles} \
  --tippecanoe ${tippecanoe} \
  --n-jobs 10 \
  --threads 10
```

| Parameter      | Required | Type   | Description                                                        |
|----------------|----------|--------|--------------------------------------------------------------------|
| `--fic-dir`    | required | string | Path to the input directory containing FICTURE2 output             |
| `--out-dir`    | required | string | Path to the output directory for PMTiles and web tiles             |
| `--id`         | required | string | Dataset ID used for naming outputs and metadata                    |
| `--makefn`     |          | string | File name for the generated Makefile (default: `run_cartload2.mk`) |
| `--spatula`    |          | string | Path to the `spatula` binary (default: `spatula`)                  |
| `--pmtiles`    |          | string | Path to the `pmtiles` binary (default: `pmtiles`)                  |
| `--tippecanoe` |          | string | Path to the `tippecanoe` binary (default: `tippecanoe`)            |
| `--n-jobs`     |          | int    | Number of parallel jobs (default: `1`)                             |
| `--threads`    |          | int    | Number of threads per job (default: `1`)                           |