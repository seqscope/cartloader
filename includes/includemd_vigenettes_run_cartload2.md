Package the SGE and the `FICTURE` results into PMTiles and a `catalog.yaml` for CartoScope. See more details in its [Reference page](../docs/reference/run_cartload2.md).

```bash
cartloader run_cartload2 \
  --makefn run_cartload2.mk \
  --fic-dir ./ficture2 \
  --out-dir ./cartload2 \
  --id ${DATA_ID} \
  --n-jobs ${n_jobs} \
  --threads ${n_jobs}
```

!!! tip "Packaging without factor analysis"
    To host only the SGE (no factors), use `cartloader run_together ... --no-ficture`, which tiles the transcripts and packages them without training any model.

<!--parameter-start-->
| Parameter      | Required | Type   | Description |
|----------------|----------|--------|-------------|
| `--fic-dir`    | required | string | `FICTURE` output directory from `run_ficture2`; must contain `ficture.params.json` (the SGE paths are read from it) |
| `--out-dir`    | required | string | Output directory for the PMTiles, assets, and `catalog.yaml` |
| `--id`         | required | string | Dataset ID used in the catalog and output metadata (no whitespace) |
| `--makefn`     |          | string | File name for the generated Makefile (default: `run_cartload2.mk`) |
| `--n-jobs`     |          | int    | Number of parallel jobs (default: `1`) |
| `--threads`    |          | int    | Number of threads per job (default: `4`) |
| `--spatula`, `--tippecanoe`, `--pmtiles` | | string | Tool paths (default: the copies built under `submodules/`) |
<!--parameter-end-->
