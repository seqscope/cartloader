# My platform is not in the supported list. Can I still use `CartLoader`?

Yes.  
If your data platform is not one of the predefined platform options, you can still use `CartLoader` by converting your input into unified SGE format with:

```bash
cartloader sge_convert --platform generic ...
```

## How to apply `CartLoader` with unsupported platforms

1. Prepare a transcript-level CSV/TSV with core fields:
   - X coordinate
   - Y coordinate
   - feature name (typically gene)
   - count (or one row per transcript)
2. Run `sge_convert` in generic mode.
3. Set CSV column names as needed using `--csv-colname-x`, `--csv-colname-y`, `--csv-colname-feature-name`, and `--csv-colnames-count`.
4. Set spatial scaling using `--units-per-um` (or your equivalent known scale).
5. Continue with standard downstream modules (`run_ficture2`, `run_cartload2`, uploads).

## Example

```bash
cartloader sge_convert \
  --platform generic \
  --in-csv /path/to/input.tsv \
  --csv-delim $'\t' \
  --csv-colname-x x \
  --csv-colname-y y \
  --csv-colname-feature-name gene \
  --csv-colnames-count count \
  --units-per-um 1 \
  --out-dir /path/to/out/sge
```

## Notes

- Start with a small subset first to validate formatting and scaling.
- If your file has comments, delimiters, or headers that differ from defaults, set the corresponding `--csv-*` options explicitly.

See also:

- [SGE Format Conversion reference](../reference/sge_convert.md)
- [How do I define scaling (`--units-per-um`) for `sge_convert`?](./define_units_per_um.md)
