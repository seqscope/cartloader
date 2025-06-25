Convert the raw input to the unified SGE format. See more details in [Reference page](../docs/reference/sge_convert.md).

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform ${PLATFORM} \
  --in-mex ./raw \
  --in-parquet ./raw/tissue_positions.parquet \
  --scale-json ./raw/scalefactors_json.json \
  --out-dir ./sge \
  --exclude-feature-regex '^(BLANK|NegCon|NegPrb)' \
  --sge-visual \
  --spatula ${spatula} \
  --n-jobs 10
```

| Parameter                 | Required | Type   | Description                                                                                                                                             |
|---------------------------|----------|--------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--makefn`                |          | string | File name for the generated Makefile                                                                                                                    |
| `--platform`              | required | string | Platform (options: "10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "generic")|
| `--in-mex`                | required | string | Path to the input MEX directory containing gene × barcode matrix                                                                                        |
| `--in-parquet`            | required | string | Path to the `tissue_positions.parquet` file with spatial barcode metadata                                                                               |
| `--scale-json`            |          | string | Path to the `scalefactors_json.json` file with coordinate scaling information. Alternatively, user could define the scale directly use `--units-per-um` |
| `--out-dir`               | required | string | Output directory for the converted SGE files                                                                                                            |
| `--exclude-feature-regex` |          | regex  | Pattern to exclude control features (e.g., BLANK, NegCon, NegPrb)                                                                                       |
| `--sge-visual`            |          | flag   | Enable SGE visualization step (generates diagnostic image)                                                                                              |
| `--spatula`               |          | string | Path to the spatula binary                                                                                                                              |
| `--n-jobs`                |          | int    | Number of parallel jobs for processing                                                                                                                  |

