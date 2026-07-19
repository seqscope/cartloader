# Platform: Illumina StrataMap (`--platform illumina`)

Built-in profile for **Illumina StrataMap**. Ingest runs `sge_convert` (the `illumina` preset). Like SeqScope, StrataMap is a sequencing-based spatial platform whose input is a **MEX matrix** (`barcodes.tsv.gz` / `features.tsv.gz` / `matrix.mtx.gz`); the spatial coordinates are **embedded in the barcode string**, so no separate positions/scale-factor file is needed. Individual per-molecule transcripts are not available — the transcript TSV is derived from the MEX.

## Expected inputs

| Role | What it is | Column / config | Required |
|------|-----------|-----------------|:---:|
| MEX matrix | The StrataMap MEX **directory** (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`) | `mex` (alias `mex_dir`) | ✅ |
| Cell boundaries | Cell boundary polygons (CSV) → enables **cell-level** analysis | `boundaries` (aliases `cell_boundary`, `cell_boundaries`) / `--in-cell-boundary` | optional |
| Cell centroids (xy) | Per-cell centroids → adds the spatial cluster scatter | `xy` (alias `cell_xy`) / `--in-cell-xy` | optional |

The `mex` role must be a **single directory** (StrataMap ships the three files together) — it is passed to `sge_convert --in-mex`. `sge_convert` reads the coordinates out of the barcode with these StrataMap defaults: barcode delimiter `:`, X = field 3, Y = field 2, `units-per-um = 1000` (barcodes encode nanometers). Override them per run through `ingest.extra_flags` in a `--config`; see [`sge_convert`](../sge_convert.md).

Default boundary/centroid column names (overridable — see below): boundaries `cell_id` / `vertex_x` / `vertex_y`; xy `cell_id` / `X` / `Y`.

**Defaults:** FICTURE `width=12`, `n_factor=24,48,96`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`.

!!! note "How the MEX is supplied"
    The StrataMap MEX is provided as the **`mex` role** (sample sheet / config), *not* via `--in-dir`
    auto-detection — StrataMap has no standard directory layout, so this profile does not assume where
    the MEX sits. There is no single-sample `--in-mex` shortcut; use a one-row sample sheet or a `--config`.

---
## Cell-level analysis (boundaries)

StrataMap has no per-transcript cell assignment, so — exactly like MERSCOPE — supplying **cell boundaries** turns on cell-level decode: after `sge_convert` builds the `X, Y, gene, count` transcript, ingest runs `spatula tsv-add-cell-id` to assign each transcript to the polygon that contains it (point-in-polygon; transcripts inside no cell → `UNASSIGNED`), appending a `cell_id` column. The `cartloader` cell factor is then decoded via `run_ficture2_multi_cells`, and the boundaries render cell polygons. `xy` (optional) adds the per-cell spatial cluster scatter.

Without boundaries the run stays **pixel-level** (FICTURE → packaging).

!!! warning "Coordinate units"
    The boundary polygon coordinates must be in **microns**, matching the barcode-derived transcript
    coordinates — `tsv-add-cell-id` does not rescale them.

---
## 1. Single sample

Give a one-row sample sheet (or a `--config`) pointing `mex` at the matrix directory; add `boundaries` for cell-level analysis:

```
id        mex                         boundaries
MYSAMPLE  /data/stratamap/MYSAMPLE    /data/stratamap/MYSAMPLE/bounds.csv
```

```bash
cartloader run_together --platform illumina \
    --samples samples.tsv --out-dir OUT --width 12 --n-factor 48
```

Boundary/centroid column names can be overridden per run with `--colname-boundary-cell/-x/-y` and `--colname-xy-cell/-x/-y`.

---
## 2. Multi-sample (sample sheet)

One joint FICTURE model across samples sharing `--out-dir`; a mix of boundary and pixel-only samples is fine (cell analysis runs for the samples that supply `boundaries`):

```
id      mex                          boundaries                       xy
s1      /data/stratamap/s1           /data/stratamap/s1/bounds.csv    /data/stratamap/s1/cells.csv
s2      /data/stratamap/s2
```

```bash
cartloader run_together --platform illumina \
    --samples samples.tsv --out-dir OUT --width 12 --n-factor 48 -j 8
```

---
## 3. Full config (JSON)

```jsonc
// run.json
{
  "platform": "illumina", "out_dir": "OUT",
  "resources": { "n_jobs": 8, "threads": 16 },
  "samples": [
    { "id": "s1", "mex": "/data/stratamap/s1", "boundaries": "/data/stratamap/s1/bounds.csv" },
    { "id": "s2", "mex": "/data/stratamap/s2" }
  ],
  "ficture": [ { "id": "denovo", "mode": "train", "width": "12", "n_factor": "24,48,96" } ],
  // override barcode parsing if a dataset differs from the StrataMap defaults:
  "ingest": { "extra_flags": ["--sge-visual", "--bcd-delim :", "--icol-bcd-x 3", "--icol-bcd-y 2", "--units-per-um 1000"] }
}
```
```bash
cartloader run_together --config run.json
```

---
## Notes

- **Coordinates come from the barcode.** No `--pos-parquet` / `--scale-json` is used; if a run's barcode format differs, override the parsing flags via `ingest.extra_flags`.
- **Images:** StrataMap has no standard morphology image; attach any histology you have as a generic image (see [Image Modalities](../run_together_images.md)).

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
- [`sge_convert`](../sge_convert.md) — barcode/MEX parsing details.
- Related sequencing-based platform: [SeqScope starter tutorial](../../vignettes/subregion_tutorials/seqscope.md).
