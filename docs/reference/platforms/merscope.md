# Platform: MERSCOPE / MERFISH (`--platform merfish`)

Built-in profile for Vizgen MERSCOPE / MERFISH. Ingest runs `sge_convert` (the `vizgen_merscope` preset). Inputs are **individual files that may be arbitrarily named and need not share a directory**, so you can point at each one explicitly — or point `--in-dir` at a standard export directory and let the profile auto-detect them.

## Expected inputs

| Role | Standard file (auto-detected in `--in-dir`) | CLI flag / column | Expected columns (defaults) |
|------|---------------------------------------------|-------------------|-----------------------------|
| Transcripts (raw) | `detected_transcripts.csv[.gz]` | `--in-transcript` / `raw_transcript` | `global_x`, `global_y`, `gene` (`global_z` optional) |
| Cell centroids (xy) | `cell_metadata.csv` | `--in-cell-xy` / `xy` | `center_x`, `center_y`; the unnamed first column holds the cell id |
| Cell boundaries | `cell_boundaries.csv` | `--in-cell-boundary` / `boundaries` | `cell_id`, `vertex_x`, `vertex_y` |
| Cell×gene matrix | `cell_by_gene.csv` | `--in-cellxgene` / `cellxgene` | first column = cell id; remaining columns = genes (`Blank*` dropped) |

**Column names may be overridden per file** when they differ: `--colname-transcript-x/-y/-feature/-count`, `--colname-xy-cell/-x/-y` (use `--colname-xy-cell ''` for the unnamed index column), and `--colname-boundary-cell/-x/-y`.

**Defaults:** FICTURE `width=12`, `n_factor=12,24,48`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`. Coordinates are already in microns and are not shifted, so the files stay aligned.

---
## 1. Single sample (direct CLI)

Point at each file explicitly:

```bash
cartloader run_together --platform merfish \
    --in-transcript    /data/detected_transcripts.csv \
    --in-cell-xy       /data/cell_metadata.csv \
    --in-cell-boundary /data/cell_boundaries.csv \
    --id MYSAMPLE --out-dir OUT --width 12 --n-factor 24
```

Or point at a standard export directory:

```bash
cartloader run_together --platform merfish --in-dir /data/MYSAMPLE --out-dir OUT
```

`--in-dir` auto-detects `detected_transcripts.csv[.gz]`, `cell_metadata.csv`, `cell_boundaries.csv`, and `cell_by_gene.csv`.

Only a transcript source is required; with just it, the run is **pixel-level** (FICTURE → packaging). Adding a cell-count source (below) enables **cell-level** decode.

---
## 2. Multi-sample (sample sheet)

For standard directories, one `in_dir` per row; a per-row role column overrides an auto-detected file:

```
id      in_dir              cell_boundary                cellxgene
s1      /data/s1
s2      /data/s2            /data/s2/bounds.csv
s3                                                       /data/s3/cell_by_gene.csv
```

```bash
cartloader run_together --platform merfish --samples samples.tsv --out-dir OUT \
    --width 12 --n-factor 24 -j 8
```

(For a sample with no `in_dir`, give `raw_transcript` + the role columns explicitly.)

---
## 3. Full config (JSON)

```jsonc
// run.json
{
  "platform": "merfish", "out_dir": "OUT",
  "resources": { "n_jobs": 8, "threads": 16 },
  "samples": [
    { "id": "s1", "in_dir": "/data/s1" },
    { "id": "s2", "raw_transcript": "/data/s2/tx.csv", "cellxgene": "/data/s2/cell_by_gene.csv", "xy": "/data/s2/meta.csv" }
  ],
  "ficture": [ { "id": "denovo", "mode": "train", "width": "12", "n_factor": "12,24,48" } ]
}
```

---
## Cell-level analysis: three sources

MERSCOPE transcripts carry **no per-transcript cell assignment**, so cell identity comes from one of three sources. Cell analysis runs when **any** is present:

| Source | How it works | Column / flag |
|--------|--------------|---------------|
| **Cell boundaries** | After `sge_convert`, ingest runs `spatula tsv-add-cell-id` to assign each transcript to the polygon that contains it (point-in-polygon; transcripts inside no cell → `UNASSIGNED`), appending a `cell_id` column. Cells are then clustered from those per-transcript counts. | `boundaries` |
| **Cell×gene matrix** | `cell_by_gene.csv` is converted to a MEX directory (`cartloader convert_cellxgene`); cells are clustered directly from that matrix. Works **without** boundaries. | `cellxgene` |
| **Existing `cell_id` column** (rare) | If the transcript CSV already carries a cell-id column, name it with `--colname-transcript-cell <name>`; it is carried to transcript column 5 and used directly (no `tsv-add-cell-id`). | `--colname-transcript-cell` |

`--in-cell-xy` / `xy` adds the per-cell spatial cluster scatter (and cell coordinates for MEX-based clustering).

### Source precedence

When more than one source is present, the **clustering** source order is:

**cellxgene MEX → transcript `cell_id` → boundary-derived**

- A `cellxgene` MEX, when present, always drives clustering (even alongside boundaries).
- Boundaries, when present, are **still used to render cell polygons**, regardless of which source drives clustering.
- `--colname-transcript-cell` skips `tsv-add-cell-id` but does **not** override a cellxgene MEX.

### Mixed joint runs

A joint run may mix sources across samples — e.g. some samples with boundaries, some with only a `cell_by_gene.csv`. Each sample's cell-count source is resolved **per sample**, and all are decoded against the shared model and packaged from one call:

```
id      in_dir              cell_boundary          cellxgene
sA      /data/sA            /data/sA/bounds.csv
sB      /data/sB                                   /data/sB/cell_by_gene.csv
sC      /data/sC            /data/sC/bounds.csv    /data/sC/cell_by_gene.csv
```

Here `sA` clusters from boundaries, `sB` and `sC` from their MEX (with `sC`'s boundaries still rendering polygons).

!!! tip "Nearest-cell fallback"
    Assign an outside-all-cells transcript to the closest cell within a distance with
    `{ "ingest": { "assign_cell_id": { "expand_um": <µm> } } }` (boundary-based decode only).

---
## Images

MERSCOPE morphology images (DAPI, protein, poly-T, …) have non-standard paths and need the mosaic transform, so declare them explicitly. The MERSCOPE transform column is `merfish_csv` → `import_image --micron2pixel-csv`; the profile applies `--shrink-factor 5.0 --high-memory` for the large mosaics.

Image sheet (`--images images.tsv`):

```
sample    type      source            merfish_csv           color
MYSAMPLE  dapi      /p/dapi.tif       /p/transform.csv      -
MYSAMPLE  protein   /p/protein.png    /p/transform.csv      008A00
```

Or inline for a single sample:

```bash
cartloader run_together --platform merfish \
  --in-transcript /p/molecules.csv --in-cell-boundary /p/polys.csv \
  --image type=dapi,source=/p/dapi.tif,merfish_csv=/p/transform.csv \
  --image type=protein,source=/p/protein.png,color=008A00,transform=/p/transform.csv \
  --id MYSAMPLE --out-dir OUT
```

See [Image Modalities](../run_together_images.md) for the full field set and type registry (`dapi`, `protein`, `polyt`, `cellbound1..3`, `hne`, …).

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
- Starter tutorial: [MERSCOPE](../../vignettes/subregion_tutorials/merscope.md)
