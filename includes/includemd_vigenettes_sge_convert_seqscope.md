Convert the raw input to the unified SGE format. See more details in [SGE Format Conversion](../docs/reference/sge_convert.md).

<!-- ```bash
cartloader sge_convert \
    --makefn sge_convert.mk \           # (optional) file name of the output make file
    --platform generic \                # use the 'generic' platform parser (adapt as needed for others like 10x_visium_hd, seqscope etc.)
    --in-csv ./input.tsv.gz \           # path to the input.tsv.gz containing raw transcript-indexed SGE
    --csv-colnames-count gn \           # column name for expression counts in the input file (use 'gn' for unique counts in the example data)
    --csv-colname-feature-name gene \   # column name for gene symbols in the input file
    --units-per-um 1000.0 \             # scale to convert coordinates to microns (the example input data is in nanometers, use 1000.0 since 1000 nm = 1 µm)
    --out-dir ./sge \                   # path to output directory where the unified SGE will be saved
    --colnames-count count  \           # output column name for expression count
    --sge-visual \                      # (optional) enable SGE visualization step
    --spatula ${spatula} \              # (optional) path to the spatula binary
    --n-jobs 10                         # (optional) number of parallel jobs for processing
``` -->

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform ${PLATFORM} \
  --in-mex ./raw \
  --units-per-um ${SCALE} \
  --out-dir ./sge \
  --exclude-feature-regex '^(BLANK|NegCon|NegPrb)' \
  --sge-visual \
  --spatula ${spatula} \
  --n-jobs 10
```

| Parameter                 | Required | Type   | Description                                                                                                                                       |
|---------------------------|----------|--------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| `--platform`              | required | string | Platform (options: "10x_visium_hd", "seqscope", "10x_xenium", "bgi_stereoseq", "cosmx_smi", "vizgen_merscope", "pixel_seq", "nova_st", "generic") |
| `--in-mex`                | required | string | Path to the input MEX directory containing gene × barcode matrix                                                                                  |
| `--units-per-um`          | required | float  | Scale to convert coordinates to microns (default: `1.0`)                                                                                          |
| `--out-dir`               | required | string | Output directory for the converted SGE files                                                                                                      |
| `--makefn`                |          | string | File name for the generated Makefile (default: `sge_convert.mk`)                                                                                  |
| `--exclude-feature-regex` |          | regex  | Pattern to exclude control features                                                                                                               |
| `--sge-visual`            |          | flag   | Enable SGE visualization step (generates diagnostic image) (default: `FALSE`)                                                                     |
| `--spatula`               |          | string | Path to the spatula binary (default: `spatula`)                                                                                                   |
| `--n-jobs`                |          | int    | Number of parallel jobs for processing (default: `1`)                                                                                             |