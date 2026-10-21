# Platform: 10x Visium HD (`--platform 10x_visium_hd`)

Built-in profile. Ingest runs `sge_convert` on the 2 µm bins. Point `--in-dir` at the Space Ranger **`outs/`** directory. Bin positions are converted to µm with `microns_per_pixel` from the 2 µm `scalefactors_json.json`; the profile's scale settings only set the resolution of the rendered raster images (2 µm per pixel), not the coordinates.

## Expected inputs (under `--in-dir`)

| Purpose | Path under `--in-dir` | Used by |
|---------|-----------------------|---------|
| Scale factors | `binned_outputs/square_002um/spatial/scalefactors_json.json` | ingest (`--scale-json`); H&E scale; centroid scaling |
| Count matrix | `binned_outputs/square_002um/filtered_feature_bc_matrix/` | ingest (`--in-mex`) |
| Bin positions | `binned_outputs/square_002um/spatial/tissue_positions.parquet` | ingest (`--pos-parquet`) |
| Square layers | `binned_outputs/square_008um`, `binned_outputs/square_016um` | `spaceranger-sq008`, `spaceranger-sq016` layers (imported if present) |
| Segmented cells | `segmented_outputs/` | `spaceranger` cell layer (imported if present) |
| Segmented cell matrix | `segmented_outputs/filtered_feature_cell_matrix` (`mex` role) | `cartloader` cell analysis |
| Segmented cell polygons | `segmented_outputs/cell_segmentations.geojson` (`boundaries` role) | `cartloader` cell analysis (optional) |
| H&E image | not auto-detected: the sample's `hne` column / config key | `hne` layer |

**Defaults:** FICTURE `width=12`, `n_factor=24,48,96`, single-molecule mode, `decode_scale=2`; packaging with `--use-pmpoint --bin-count 500 --sge-scale 2`; `sge_convert` with `--sge-visual --scale-xy 2.0`. These three scale settings render the ingest preview (`xy.png`), the FICTURE pixel-decode images, and the packaged SGE raster at 2 µm per pixel.

---
## Cell layers

Two cell-level layers can come out of a Space Ranger run that includes segmentation:

| Layer | What it is | Runs when |
|-------|-----------|-----------|
| `spaceranger` | Space Ranger's own segmented cells and clusters, imported by [`import_visiumhd_cell`](../import_cell.md) in the images stage | `segmented_outputs/` exists |
| `cartloader` | Cells clustered from the segmented cell matrix and decoded against the FICTURE model by [`run_ficture2_multi_cells`](../run_ficture2_multi_cells.md) | the sample has the segmented cell matrix (`mex` role); on a joint run, **every** sample must have it (`require_all`) |

For the `cartloader` analysis, cell centroids are derived from `cell_segmentations.geojson` when it is present: the polygon ids (property `cell_id`) are formatted as `cellid_{:09d}-1` to match the matrix barcodes and rescaled to µm with `microns_per_pixel` from the scale factors JSON. Its polygons are not rendered by this analysis (the `spaceranger` layer renders them).

The square layers (`spaceranger-sq008`, `spaceranger-sq016`) are imported by [`import_visiumhd_square`](../import_square.md) when those directories exist.

---
## 1. Single sample (direct CLI)

```bash
cartloader run_together --platform 10x_visium_hd \
    --in-dir /data/visiumhd/sample/outs --out-dir OUT \
    --width 12 --n-factor 48 -j 8 --threads 8
```

Square (8/16 µm) and cell layers are added automatically when present. The bins sit on a 2 µm lattice; `--jitter-xy 0.8` spreads each bin's molecules across the bin (see [Coordinate jitter](../run_together_inputs.md#coordinate-jitter-at-ingest)).

**H&E** is not inside the standard bin layout, so name it in the sample's `hne` field — for a single sample, a one-row sample sheet:

```
in_dir                         hne
/data/visiumhd/sample/outs     /data/visiumhd/sample/he.tif
```

```bash
cartloader run_together --platform 10x_visium_hd \
    --samples sample.tsv --out-dir OUT --width 12 --n-factor 48
```

Without an `id` column the lone sample is `rep1`, exactly as with `--in-dir`. Give the `hne` path as absolute (it is not resolved against `in_dir`).

!!! warning "`--image type=hne` does not use the scale factors"
    An H&E given through `--image` or `--images` is an ordinary `rgb` image: it does not read `scalefactors_json.json`, so it needs its own bounds (`georef_plain=true,um_per_pixel=<µm per pixel>`), or the image import stops with an error. Only the sample's `hne` field gets the profile's H&E handling.

---
## 2. Multi-sample (sample sheet)

One joint model; add a per-sample `hne` column for the H&E layer:

```bash
cartloader run_together --platform 10x_visium_hd \
    --samples samples.tsv --out-dir OUT --width 12 --n-factor 48 -j 8
```

```
id      in_dir                       hne
s1      /data/visiumhd/s1/outs       /data/visiumhd/s1/he.tif
s2      /data/visiumhd/s2/outs       /data/visiumhd/s2/he.tif
```

---
## 3. Full config (JSON)

```bash
cartloader run_together --platform 10x_visium_hd --config run.json
```
```jsonc
// run.json
{
  "platform": "10x_visium_hd", "out_dir": "OUT",
  "resources": { "n_jobs": 8, "threads": 16 },
  "samples": [
    { "id": "s1", "in_dir": "/data/visiumhd/s1/outs", "hne": "/data/visiumhd/s1/he.tif" },
    { "id": "s2", "in_dir": "/data/visiumhd/s2/outs", "hne": "/data/visiumhd/s2/he.tif" }
  ],
  "ficture": [ { "id": "denovo", "mode": "train", "width": "12", "n_factor": "24,48,96" } ]
}
```

---
## Notes

- **H&E** is an `rgb` (passthrough) image imported as the `hne` layer with `--georef-plain`: its upper-left corner is placed at (0,0) and its scale is `microns_per_pixel` from the 2 µm `scalefactors_json.json` (read with `jq` at run time). See [Image Modalities](../run_together_images.md#georeferencing).
- Square (8/16 µm) and `spaceranger` cell imports are wired by the profile's `squares` / `cell_import` blocks and appended to the catalog when the directories exist.

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
- Step-by-step (manual modules): [Visium HD end-to-end tutorial](../../vignettes/pipelines/visiumhd.md)
