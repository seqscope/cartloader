
Convert the raw input to the unified SGE format. See more details in its [Reference page](../docs/reference/sge_convert.md).

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform ${PLATFORM} \
  --in-csv ./input.tsv.gz \
  --units-per-um ${SCALE} \
  --out-dir ./sge \
  --exclude-feature-regex '^(BLANK|Neg|Intergenic|Deprecated|Unassigned)' \
  --sge-visual \
  --n-jobs ${n_jobs}
```

<!--parameter-start-->
| Parameter                    | Required | Type   | Description |
|------------------------------|----------|--------|-------------|
| `--platform`                 | required | string | Input platform, which sets the default input columns and delimiter (options: `10x_visium_hd`, `seqscope`, `illumina`, `10x_xenium`, `bgi_stereoseq`, `cosmx_smi`, `vizgen_merscope`, `pixel_seq`, `nova_st`, `generic`) |
| `--in-csv`                   | required | string | Path to the input transcript CSV/TSV file |
| `--out-dir`                  | required | string | Output directory for the converted SGE files |
| `--units-per-um`             |          | float  | Coordinate units per µm in the input (default: `1`, i.e. already in µm; `1000` for the Seq-Scope raw TSV; detected from the barcodes for `illumina`) |
| `--csv-colname-feature-name` |          | string | Gene-name column in the input (default: set by `--platform`) |
| `--csv-colname-count`        |          | string | Count column in the input (default: set by `--platform`; without one, each row counts 1) |
| `--makefn`                   |          | string | File name for the generated Makefile (default: `sge_convert.mk`) |
| `--exclude-feature-regex`    |          | regex  | Pattern of feature names to drop (e.g. control probes) |
| `--sge-visual`               |          | flag   | Draw the converted SGE as an image (`xy.png`) for a visual check |
| `--n-jobs`                   |          | int    | Number of parallel jobs (default: `1`) |
<!--parameter-end-->
