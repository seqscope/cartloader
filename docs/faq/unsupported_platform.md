# My platform is not in the supported list. Can I still use `CartLoader`?

Yes. Any platform that can produce a transcript-level table works.

## 1. Prepare a transcript table

A transcript-level CSV/TSV with these fields:

- X coordinate
- Y coordinate
- feature name (typically gene)
- count (or one row per transcript)

## 2. Run `run_together` with the generic profile (recommended)

```bash
cartloader run_together \
  --platform generic \
  --raw-transcript /path/to/input.tsv.gz \
  --colname-transcript-x x \
  --colname-transcript-y y \
  --colname-transcript-feature gene \
  --colname-transcript-count count \
  --units-per-um 1 \
  --out-dir /path/to/out/my_sample \
  --width 18 --n-factor 12
```

The [Generic tutorial](../vignettes/subregion_tutorials/generic.md) walks through this on real data, including a comma-separated file and a cell-ID column. `--raw-transcript` sends the table through SGE conversion; the `--colname-transcript-*` options name its columns, and `--units-per-um` sets how many coordinate units make one micrometer. See the [Generic platform page](../reference/platforms/generic.md) for every option, including a cell-ID column for cell analyses.

## Alternative: convert with `sge_convert` yourself

To control the conversion directly, run [`sge_convert`](../reference/sge_convert.md) in generic mode, then continue with the [step-by-step](../vignettes/step_by_step.md) commands (`run_ficture2`, `run_cartload2`, uploads):

```bash
cartloader sge_convert \
  --platform generic \
  --in-csv /path/to/input.tsv \
  --csv-delim $'\t' \
  --csv-colname-x x \
  --csv-colname-y y \
  --csv-colname-feature-name gene \
  --csv-colname-count count \
  --units-per-um 1 \
  --out-dir /path/to/out/sge
```

## Notes

- Start with a small subset first to validate formatting and scaling.
- If your file has comments, delimiters, or headers that differ from defaults, set the corresponding options explicitly.

See also:

- [How do I define scaling (`--units-per-um`)?](./define_units_per_um.md)
- [SGE Format Conversion reference](../reference/sge_convert.md)
