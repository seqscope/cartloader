Convert the raw input to the unified SGE format. See more details in [SGE Format Conversion](../docs/reference/sge_convert.md).

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform ${PLATFORM} \
  --in-mex ./raw \
  --units-per-um ${SCALE} \
  --out-dir ./sge \
  --exclude-feature-regex '^(BLANK|NegCon|NegPrb)' \
  --sge-visual \
  --n-jobs ${n_jobs}
```

| Parameter                 | Required | Type   | Description |
|---------------------------|----------|--------|-------------|
| `--platform`              | required | string | Input platform (options: `10x_visium_hd`, `seqscope`, `illumina`, `10x_xenium`, `bgi_stereoseq`, `cosmx_smi`, `vizgen_merscope`, `pixel_seq`, `nova_st`, `generic`) |
| `--in-mex`                | required | string | Input MEX directory holding `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` |
| `--out-dir`               | required | string | Output directory for the converted SGE files |
| `--units-per-um`          |          | float  | Coordinate units per µm in the input (default: `1`). Seq-Scope MEX coordinates are in nanometers, so pass `1000` |
| `--icols-mtx`             |          | int    | Which of the matrix's count columns to use, 1-based (default: `1`, the "Gene" count of Seq-Scope data) |
| `--makefn`                |          | string | File name for the generated Makefile (default: `sge_convert.mk`) |
| `--exclude-feature-regex` |          | regex  | Pattern of feature names to drop (e.g. control probes) |
| `--sge-visual`            |          | flag   | Draw the converted SGE as an image (`xy.png`) for a visual check |
| `--n-jobs`                |          | int    | Number of parallel jobs (default: `1`) |
