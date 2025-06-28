
Convert the raw input to the unified SGE format. See more details in [Reference page](../docs/reference/sge_convert.md).

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform ${PLATFORM} \
  --in-csv ./input.tsv.gz \
  --units-per-um ${SCALE} \
  --out-dir ./sge \
  --exclude-feature-regex '^(BLANK|NegCon|NegPrb)' \
  --sge-visual \
  --spatula ${spatula} \
  --n-jobs 10
```

| Parameter                 | Required | Type   | Description                                                                                                                                       | 
|---------------------------|----------|--------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| `--platform`              | required | string | Platform (options: "`10x_visium_hd`", "`seqscope`", "`10x_xenium`", "`bgi_stereoseq`", "`cosmx_smi`", "`vizgen_merscope`", "`pixel_seq`", "`generic`") | 
| `--in-csv`                | required | string | Path to the input TSV/CSV file                                     |
| `--units-per-um`          | required | float  | Scale to convert coordinates to microns (default: `1.0`)                                                                                           | 
| `--out-dir`               | required | string | Output directory for the converted SGE files                                                                                                      | 
| `--makefn`                |          | string | File name for the generated Makefile (default: `sge_convert.mk`)                                                                                  |
| `--exclude-feature-regex` |          | regex  | Pattern to exclude control features                                                                                                               |
| `--sge-visual`            |          | flag   | Enable SGE visualization step (generates diagnostic image) (default: `FALSE`)                                                                     |
| `--spatula`               |          | string | Path to the spatula binary (default: `spatula`)                                                                                                   |
| `--n-jobs`                |          | int    | Number of parallel jobs for processing (default: `1`)                                                                                             |