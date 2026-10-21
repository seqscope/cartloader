Package the SGE and the `FICTURE` results into PMTiles and a `catalog.yaml` for CartoScope. See more details in [Reference page](../docs/reference/run_cartload2.md).

Use the same mount as the `FICTURE` step, so that the paths recorded in `/out/ficture2/ficture.params.json` resolve.

```bash
docker run --rm \
  -v $(pwd)/out:/out \
  hyunminkang/cartloader:latest \
  run_cartload2 \
    --makefn run_cartload2.mk \
    --fic-dir /out/ficture2 \
    --out-dir /out/cartload2 \
    --id ${DATA_ID} \
    --n-jobs ${n_jobs} \
    --threads ${n_jobs}
```
