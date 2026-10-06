# Platform: CosMx SMI (`--platform cosmx_smi`)

Built-in profile. Unlike the others, ingest does **not** run `sge_convert`: a dedicated [`reformat_cosmx`](../reformat_cosmx.md) step reads the three raw CSVs and, converting global-pixel coordinates to microns (`0.12028 µm/px`), writes the transcript TSV plus the cell-metadata (xy) and polygon (boundaries) files that FICTURE and cell decode consume. Point `--in-dir` at the AtoMx / CosMx SMI flat-file export directory; filenames are matched by glob (exactly one file must match each).

## Expected inputs (under `--in-dir`)

| Purpose | Default glob | reformat flag | Produces |
|---------|--------------|:---:|----------|
| Transcripts | `*tx_file.csv.gz` (cols `fov`, `cell_ID`, `target`, `x_global_px`, `y_global_px`, `z`) | `--tx` | `*.transcripts.tsv.gz` |
| Cell centroids | `*metadata*.csv.gz` (cols `fov`, `cell_ID`, `CenterX_global_px`, `CenterY_global_px`) | `--meta` | `*.metadata.csv.gz` (xy, cols `X`/`Y`) |
| Cell boundaries | `*polygons.csv.gz` (cols `fov`, `cellID`, `x_global_px`, `y_global_px`) | `--poly` | `*.polygons.csv.gz` |

Cell ids are formed as `<fov>_<cell_ID>`; transcripts with `cell_ID = 0` are kept but tagged `UNASSIGNED`, and `System*` control probes are dropped at ingest. The `cartloader` cell factor is decoded from the metadata + polygon files via `run_ficture2_multi_cells`.

**Defaults:** FICTURE `width=12`, `n_factor=12,24,48`, `min_ct_per_unit_hexagon=100`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`.

---
## 1. Single sample (direct CLI)

```bash
cartloader run_together --platform cosmx_smi \
    --in-dir /data/cosmx/export --out-dir OUT \
    --width 12 --n-factor 24 -j 8 --threads 8
```

### Overriding the input file patterns

CosMx exports don't always use the standard suffixes. Any glob can be overridden with a small JSON whose `ingest.inputs` block **deep-merges** over the profile — restate only the patterns that differ:

```json
{ "platform": "cosmx_smi", "ingest": { "inputs": { "--tx": "*tx_unique.csv.gz" } } }
```
```bash
cartloader run_together --config cosmx.json --in-dir /data/cosmx/export --out-dir OUT
```

For an export with `..._tx.csv.gz`, `..._tx_unique.csv.gz`, `..._metadata.csv.gz`, and `..._polygons.csv.gz`, only the `--tx` override above is needed (the default `--meta`/`--poly` globs already match `_metadata`/`_polygons`).

---
## 2. Multi-sample (sample sheet)

For raw exports, one `in_dir` per row (each reformatted independently, then jointly modeled):

```
id      in_dir
s1      /data/cosmx/s1
s2      /data/cosmx/s2
```

For **already-reformatted** samples, point the role columns at the produced files to **skip re-ingest**:

```
id      transcript                    cell_xy                    cell_boundary
s1      /out/s1/s1.transcripts.tsv.gz /out/s1/s1.metadata.csv.gz /out/s1/s1.polygons.csv.gz
```

```bash
cartloader run_together --platform cosmx_smi --samples samples.tsv --out-dir OUT \
    --width 12 --n-factor 24 -j 8
```

---
## 3. Full config (JSON)

```jsonc
// run.json — per-sample glob override + joint samples
{
  "platform": "cosmx_smi", "out_dir": "OUT",
  "resources": { "n_jobs": 8, "threads": 16 },
  "ingest": { "inputs": { "--tx": "*tx_unique.csv.gz" } },
  "samples": [
    { "id": "s1", "in_dir": "/data/cosmx/s1" },
    { "id": "s2", "in_dir": "/data/cosmx/s2" }
  ]
}
```
```bash
cartloader run_together --config run.json
```

---
## Notes

- **Glob discipline:** each pattern must match exactly one file — zero or multiple is an error that lists what's present. Narrow the pattern via `ingest.inputs`.
- Transcripts already carry a per-molecule `cell_ID`, so no `tsv-add-cell-id` step is needed (contrast with MERSCOPE).

## See also

- [`reformat_cosmx`](../reformat_cosmx.md) · [Specifying Inputs](../run_together_inputs.md) · [Overview](../run_together.md)
- Starter tutorial: [CosMx SMI](../../vignettes/subregion_tutorials/cosmxsmi.md)
