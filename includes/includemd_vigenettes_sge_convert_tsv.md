
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

| Parameter                 | Required | Type   | Description                                                        |
|---------------------------|----------|--------|--------------------------------------------------------------------|
| `--makefn`                |          | string | File name for the generated Makefile                               |
| `--platform`              | required | string | platform (options: "10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "generic")                             |
| `--in-csv`                | required | string | Path to the input TSV/CSV file                                     |
| `--units-per-um`          | required | float  | Conversion factor: how many input coordinate units per micron      |
| `--out-dir`               | required | string | Output directory for the generated SGE files                       |
| `--exclude-feature-regex` |          | regex  | Pattern to exclude features (e.g., BLANK, NegCon)                  |
| `--sge-visual`            |          | flag   | Enable generation of diagnostic visualization of the SGE matrix    |
| `--spatula`               |          | string | Path to the spatula binary                                         |
| `--n-jobs`                |          | int    | Number of parallel jobs for processing                             |