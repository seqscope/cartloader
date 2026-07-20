# Platform: Seq-Scope (`--platform seqscope`)

Built-in profile for Seq-Scope. The input is a **MEX triple** (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`), but in the Seq-Scope dialect: the barcode file carries the spatial coordinates, and the feature file and matrix carry **five count columns per entry** instead of one. `sge_convert --platform seqscope` (via `spatula convert-sge`) already reads that dialect; this profile wires it into the end-to-end run.

Seq-Scope has **no cell segmentation**, so the run is pixel-level throughout: ingest → FICTURE → packaging (+ histology). No cell analysis stage is planned or expected.

## Expected inputs

`--in-dir` **is the MEX directory itself** — there is no standard parent layout to search:

| File in `--in-dir` | Contents |
|--------------------|----------|
| `barcodes.tsv.gz` | `barcode`, index, …, **X (col 6)**, **Y (col 7)**, counts — coordinates in **nanometers** |
| `features.tsv.gz` | `gene_id`, `gene_name`, index, counts |
| `matrix.mtx.gz` | MatrixMarket with **5 count columns**: `gn, gt, spl, unspl, ambig` |

Histology has no standard filename, so it is always given explicitly with `--image` (see [Images](#images)).

**Defaults:** ingest with `--units-per-um 1000` (nm → µm), `--icols-mtx 2` (the gene-exon count), `--sge-visual`; FICTURE `width=12`, `n_factor=24,48,96`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`.

---
## 1. Single sample

```bash
cartloader run_together --platform seqscope \
    --in-dir /data/seqscope/sample1/mex \
    --out-dir OUT \
    --image type=hne,source=/data/seqscope/sample1/HE.tif \
    --width 12 --n-factor 24
```

The sample id defaults to the `--out-dir` basename; override with `--id`.

---
## 2. Multi-sample (sample sheet)

Each row's `in_dir` is that sample's MEX directory. Images come from a separate `--images` TSV keyed by `sample`:

`samples.tsv`

```
id     in_dir
liver1 /data/seqscope/liver1/mex
liver2 /data/seqscope/liver2/mex
```

`images.tsv`

```
sample type source
liver1 hne  /data/seqscope/liver1/HE.tif
liver2 hne  /data/seqscope/liver2/HE.tif
```

```bash
cartloader run_together --platform seqscope \
    --samples samples.tsv --images images.tsv --out-dir OUT -j 8
```

---
## Coordinates and count columns

Barcode X/Y are in **nanometers**, so ingest passes `--units-per-um 1000`. Columns 6 and 7 of `barcodes.tsv.gz` are read as X and Y (the `sge_convert` defaults `--icol-bcd-x 6 --icol-bcd-y 7`); a different layout is set through the config:

```json
{ "ingest": { "extra_flags": ["--units-per-um 1000", "--icols-mtx 2", "--icol-bcd-x 5", "--icol-bcd-y 6"] } }
```

The matrix's five count columns are, in order, `gn, gt, spl, unspl, ambig`. The profile selects **column 2** (`--icols-mtx 2`), the gene-exon count, and names it `count` in the output. To quantify on a different column — e.g. total gene counts — override the flag the same way (`"--icols-mtx 1"`).

Barcodes whose features do not match the feature file are dropped by `sge_drop_mismatches` after conversion (the `sge_convert` default for this platform).

---
## Images

Seq-Scope histology is typically an **H&E TIF that has already been georeferenced** — it carries its own CRS/geotransform, aligned to the transcript frame. The profile therefore sets `image_defaults.georeferenced = true`, and the image is tiled as-is:

```bash
cartloader image_png2pmtiles --in-img HE.tif --out-prefix .../hne --geotif2mbtiles --mbtiles2pmtiles
```

No bounds are synthesized, so `--um-per-pixel` / `--georef-plain` do not apply. If a particular image is **not** already georeferenced, opt it out per image and state its scale instead:

```bash
--image type=hne,source=/data/HE_plain.tif,georeferenced=false,georef_plain=true,um_per_pixel=0.5
```

Single-channel modalities (e.g. `type=dapi`) are colorized through `import_image` on the usual path and are not covered by `georeferenced`; give them a scale with `um_per_pixel`.

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
