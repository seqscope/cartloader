# Platform: Illumina StrataMap (`--platform illumina`)

Built-in profile for **Illumina StrataMap**. Ingest runs `sge_convert` (the `illumina` preset). Like SeqScope, StrataMap is a sequencing-based spatial platform whose input is a **MEX matrix** (`barcodes.tsv.gz` / `features.tsv.gz` / `matrix.mtx.gz`); the spatial coordinates are **embedded in the barcode string**, so no separate positions/scale-factor file is needed. Individual per-molecule transcripts are not available — the transcript TSV is derived from the MEX.

## Expected inputs

| Role | What it is | Column / config | Required |
|------|-----------|-----------------|:---:|
| MEX matrix | The StrataMap MEX **directory** (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`) | `mex` (alias `mex_dir`) | ✅ |
| Cell boundaries | Cell boundary polygons (CSV) → enables **cell-level** analysis | `boundaries` (aliases `cell_boundary`, `cell_boundaries`) / `--in-cell-boundary` | optional |
| Cell centroids (xy) | Per-cell centroids → adds the spatial cluster scatter | `xy` (alias `cell_xy`) / `--in-cell-xy` | optional |

The `mex` role must be a **single directory** (StrataMap ships the three files together) — it is passed to `sge_convert --in-mex`. `sge_convert` reads the coordinates out of the barcode with these StrataMap defaults: barcode delimiter `:`, X = field 3, Y = field 2. The **coordinate units are detected from the barcode file** (see below); override any of these per run — the units with `--units-per-um`, the rest through `ingest.extra_flags` in a `--config`; see [`sge_convert`](../sge_convert.md).

### Barcode coordinate units (nanometers or microns)

StrataMap has shipped two barcode formats, and they differ by a factor of 1000:

| Barcode | Units | `--units-per-um` |
|---|---|---|
| `SBC:433503:2393851` (older) | nanometers | `1000` |
| `SBC:686.951:4668.15` (current) | microns | `1` |

Ingest reads the barcode file and picks the right one — a fractional coordinate means microns, otherwise the coordinate magnitude decides — and prints the value it used, so **both formats run correctly without a flag**. State it yourself when you want it pinned, or when a file is unusual enough that detection should not be trusted:

```bash
cartloader run_together --platform illumina --samples samples.tsv --out-dir OUT --units-per-um 1
```

`--units-per-um` always wins over the detection; if it contradicts the file, the run warns and uses your value. The config equivalent is `{"ingest": {"units_per_um": 1}}`. Getting this wrong rescales the whole sample by 1000x, so it is worth checking the reported units in the ingest log against `gzip -cd barcodes.tsv.gz | head`.

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
  // units are detected from the barcode file; pin them (and any other barcode parsing
  // that differs from the StrataMap defaults) here:
  "ingest": { "units_per_um": 1,
              "extra_flags": ["--sge-visual", "--bcd-delim :", "--icol-bcd-x 3", "--icol-bcd-y 2"] }
}
```
```bash
cartloader run_together --config run.json
```

---
## Notes

- **Coordinates come from the barcode.** No `--pos-parquet` / `--scale-json` is used. The nanometer/micron convention is detected per dataset (see [Barcode coordinate units](#barcode-coordinate-units-nanometers-or-microns)) and can be pinned with `--units-per-um`; if a run's barcode *layout* differs, override the parsing flags via `ingest.extra_flags`.
- **Images:** StrataMap has no standard morphology image; attach any histology you have as a generic image (see [Image Modalities](../run_together_images.md)).

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
- [`sge_convert`](../sge_convert.md) — barcode/MEX parsing details.
- Related sequencing-based platform: [SeqScope starter tutorial](../../vignettes/subregion_tutorials/seqscope.md).
