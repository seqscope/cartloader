# Platform: Generic (`--platform generic`)

Built-in profile (`assets/run_together_profiles/generic.json`) for spatial data from **any source** that can be expressed as a transcript table: one row per molecule (or per spatial bin and gene) with X, Y, gene, and optionally a count and a cell id. Use it for platforms without a dedicated profile, or for data you have already converted yourself.

A sample's transcript comes in one of two forms:

- **Already ingested** — a TSV in CartLoader's transcript layout; ingest is skipped.
- **Raw** — a CSV/TSV converted by `sge_convert --platform generic` at ingest.

## Expected inputs

| Role | Auto-detected under `--in-dir` | CLI flag / column | Format |
|------|--------------------------------|-------------------|--------|
| Transcript (already ingested) | `transcripts.tsv.gz` | `--in-transcript` / `transcript` (`tsv`) | Tab-separated with a header row `X  Y  gene  count`, optionally a fifth `cell_id` column; coordinates in µm. Ingest filters, `--jitter-xy` and `--units-per-um` do not apply to it. |
| Transcript (raw) | — | `--raw-transcript` / `raw_transcript` | Defaults: tab-separated, columns `X`, `Y`, `gene`, no count column (each row counts 1), coordinates in µm. Override with `--colname-transcript-x/-y/-feature/-count`, `--units-per-um`, and `ingest.extra_flags` (e.g. `["--csv-delim ,"]` for a comma-separated file). A `.parquet` is converted to `.csv.gz` first. |
| Cell centroids | `cells.csv.gz` | `--in-cell-xy` / `xy` (`cell_xy`) | Columns `X`, `Y`, and a cell id (`--colname-cell` or `--colname-xy-cell`, else `cell_id`). |
| Cell boundaries | `cell_boundaries.csv.gz` | `--in-cell-boundary` / `boundaries` | Comma-separated; the first three columns are `cell_id`, `vertex_x`, `vertex_y` (µm), vertices grouped by cell. Used to render cell polygons. |
| Cell-feature MEX | — | `mex` (`mex_dir`) or `mex_bcd` / `mex_ftr` / `mex_mtx`; or `cellxgene` (`cell_by_gene`, `--in-cellxgene`), converted to a MEX at ingest for a raw-transcript sample | Barcodes are the cell ids. |

Only a transcript is required. With `--in-dir`, a directory without `transcripts.tsv.gz` has no transcript source and the run stops; give the file with `--in-transcript` or `--raw-transcript` instead.

**Defaults:** FICTURE `width=12`, `n_factor=12,24,48`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`; no images (attach them with [`--image` / `--images`](../run_together_images.md)).

---
## Cell-level analysis

The `cartloader` cell analysis is selected **explicitly**, for all samples of a run at once:

| Cell-count source | How to turn it on |
|-------------------|-------------------|
| A `cell_id` column in the transcript (column 5) | `--colname-cell <name>` (or top-level `colname_cell` in `--config`). The column is read by position; the name is also the default cell-id column of the `xy` file. For a raw transcript, also pass `--colname-transcript-cell <raw column>` so ingest carries that column to position 5. |
| A MEX in **every** sample | Provide `mex` (or `cellxgene`) for each sample; no flag needed. |

Giving both is an error; with neither, cell analysis is skipped and the run is pixel-level (FICTURE → packaging). When present, `xy` and `boundaries` are added to the decode (cell centroids, polygons). `--colname-boundary-*` has no effect on this platform.

---
## 1. Single sample (direct CLI)

A directory holding the already-ingested files under their standard names:

```
/data/sample1/
├── transcripts.tsv.gz       # X, Y, gene, count, cell_id  (µm)
├── cells.csv.gz             # optional: cell_id, X, Y
└── cell_boundaries.csv.gz   # optional: cell_id, vertex_x, vertex_y
```

```bash
cartloader run_together --platform generic \
    --in-dir /data/sample1 --out-dir OUT --colname-cell cell_id \
    --width 12 --n-factor 24
```

A raw comma-separated file with its own column names:

```bash
cartloader run_together --platform generic \
    --raw-transcript /data/molecules.csv \
    --colname-transcript-x x_um --colname-transcript-y y_um \
    --colname-transcript-feature gene_name \
    --config generic.json --out-dir OUT --width 12 --n-factor 24
```
```json
{ "ingest": { "extra_flags": ["--csv-delim ,"] } }
```

As on every platform, a lone `--out-dir` sample is `rep1` unless `--id` is given, and its packaged directory is `OUT/cartl/<OUT basename>-rep1/`.

---
## 2. Multi-sample (sample sheet)

One joint FICTURE model across samples sharing `--out-dir`. Here every sample has a cell-feature MEX, so the cell analysis is MEX-based:

```
id    transcript                     mex                   xy
s1    /data/s1/transcripts.tsv.gz    /data/s1/cell_mex     /data/s1/cells.csv.gz
s2    /data/s2/transcripts.tsv.gz    /data/s2/cell_mex     /data/s2/cells.csv.gz
```

```bash
cartloader run_together --platform generic \
    --samples samples.tsv --out-dir OUT --width 12 --n-factor 24 -j 8
```

---
## 3. Full config (JSON)

```jsonc
// run.json
{
  "platform": "generic", "out_dir": "OUT",
  "colname_cell": "cell_id",
  "samples": [
    { "id": "s1", "in_dir": "/data/s1" },
    { "id": "s2", "transcript": "/data/s2/tx.tsv.gz", "xy": "/data/s2/cells.csv.gz" }
  ],
  "ficture": [ { "id": "denovo", "mode": "train", "width": "12", "n_factor": "12,24,48" } ]
}
```
```bash
cartloader run_together --config run.json
```

---
## Notes

- **A different `sge_convert` preset.** To read a raw file with another platform's `sge_convert` preset (e.g. `pixel_seq`, `nova_st`), set `ingest.sge_platform` in a `--config` or a `--platform-json` file; see [Supported Platforms → Custom profiles](../supported_platforms.md#custom-profiles).
- **Coordinates must be in µm** after ingest. Boundary and centroid files are never rescaled.

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
- [`sge_convert`](../sge_convert.md) — the `generic` CSV/TSV ingest options.
- [Unsupported platform?](../../faq/unsupported_platform.md)
