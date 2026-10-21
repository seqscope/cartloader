Convert the raw input to the unified SGE format. See more details in its [Reference page](../docs/reference/sge_convert.md).

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform ${PLATFORM} \
  --in-mex ./raw \
  --pos-parquet ./raw/tissue_positions.parquet \
  --scale-json ./raw/scalefactors_json.json \
  --out-dir ./sge \
  --exclude-feature-regex '^(BLANK|NegCon|NegPrb)' \
  --sge-visual \
  --n-jobs ${n_jobs}
```

| Parameter                 | Required              | Type   | Description |
|---------------------------|-----------------------|--------|-------------|
| `--platform`              | required              | string | Input platform (options: `10x_visium_hd`, `seqscope`, `illumina`, `10x_xenium`, `bgi_stereoseq`, `cosmx_smi`, `vizgen_merscope`, `pixel_seq`, `nova_st`, `generic`) |
| `--in-mex`                | required              | string | Input MEX directory holding `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` |
| `--pos-parquet`           | required              | string | Path to `tissue_positions.parquet`, which gives each barcode's position |
| `--scale-json`            | required <sup>1</sup> | string | Path to `scalefactors_json.json`, used to convert pixel positions to µm |
| `--out-dir`               | required              | string | Output directory for the converted SGE files |
| `--makefn`                |                       | string | File name for the generated Makefile (default: `sge_convert.mk`) |
| `--exclude-feature-regex` |                       | regex  | Pattern of feature names to drop (e.g. control probes) |
| `--sge-visual`            |                       | flag   | Draw the converted SGE as an image (`xy.png`) for a visual check |
| `--n-jobs`                |                       | int    | Number of parallel jobs (default: `1`) |

<sub><sup>1</sup>: `CartLoader` needs the scale either from a JSON file (`--scale-json`, which must contain `microns_per_pixel`) or as a number (`--units-per-um`, pixels per µm).</sub>
