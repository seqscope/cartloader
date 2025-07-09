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
  --n-jobs ${n_jobs}
```

| Parameter                 | Required      | Type   | Description                                                                                                                                       |
|---------------------------|---------------|--------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| `--platform`              | required      | string | Platform (options: "`10x_visium_hd`", "`seqscope`", "`10x_xenium`", "`bgi_stereoseq`", "`cosmx_smi`", "`vizgen_merscope`", "`pixel_seq`", "`generic`") |
| `--in-mex`                | required      | string | Path to the input MEX directory containing gene × barcode matrix                                                                                  |
| `--in-parquet`            | required      | string | Path to the `tissue_positions.parquet` file with spatial barcode metadata                                                                         |
| `--scale-json`            | required <sup>1</sup> | string | Path to the `scalefactors_json.json` file for coordinate scaling (or use `--units-per-um` to specify directly)   |
| `--out-dir`               | required      | string | Output directory for the converted SGE files                                                                                                      |
| `--makefn`                |               | string | File name for the generated Makefile (default: `sge_convert.mk`)                                                                                  |
| `--exclude-feature-regex` |               | regex  | Pattern to exclude control features                                                                                                               |
| `--sge-visual`            |               | flag   | Enable SGE visualization step (generates diagnostic image) (default: `FALSE`)                                                                     |
| `--spatula`               |               | string | Path to the spatula binary (default: `spatula`)                                                                                                   |
| `--n-jobs`                |               | int    | Number of parallel jobs for processing (default: `1`)                                                                                             |

<sub><sup>1</sup>: To define the scaling factor, `cartloader` requires either a JSON file (via `--scale-json`) or a direct scale value using `--units-per-um`. When using `--scale-json`, make sure the JSON file has `microns_per_pixel` information.</sub>
