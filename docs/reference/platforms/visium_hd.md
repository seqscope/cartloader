# Platform: 10x Visium HD (`--platform 10x_visium_hd`)

Built-in profile. Ingest runs `sge_convert`. Point `--in-dir` at the Space Ranger **`outs/`** directory. Coordinates are scaled by 2 automatically (`--scale-xy 2.0` at ingest, `--sge-scale 2` at packaging).

## Expected inputs (under `--in-dir`)

| Purpose | Path under `--in-dir` |
|---------|-----------------------|
| Scale factors | `binned_outputs/square_002um/spatial/scalefactors_json.json` |
| Count matrix | `binned_outputs/square_002um/filtered_feature_bc_matrix/` |
| Bin positions | `binned_outputs/square_002um/spatial/tissue_positions.parquet` |
| Square layers | `binned_outputs/square_008um`, `binned_outputs/square_016um` (imported if present) |
| Segmented cells | `segmented_outputs/` (imported if present) |
| H&E image | Provided per sample via the `hne` column / config field; µm/pixel read from the 2 µm `scalefactors_json.json` |

**Defaults:** FICTURE `width=12`, `n_factor=24,48,96`, `decode_scale=2`; packaging with `--use-pmpoint --bin-count 500 --sge-scale 2`.

---
## 1. Single sample (direct CLI)

```bash
cartloader run_together --platform 10x_visium_hd \
    --in-dir /data/visiumhd/sample/outs --out-dir OUT \
    --width 12 --n-factor 48 -j 8 --threads 8
```

Square (8/16 µm) and segmented-cell layers are imported automatically when present.

**H&E** is not inside the standard bin layout, so name it. Single sample, inline:

```bash
cartloader run_together --platform 10x_visium_hd \
    --in-dir /data/visiumhd/sample/outs --out-dir OUT \
    --image type=hne,source=/data/visiumhd/sample/he.tif --width 12 --n-factor 48
```

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

- **Automatic scaling.** The 2 µm bin coordinates are scaled ×2 at ingest and packaging — you do not set this yourself.
- **H&E** is an `rgb` (passthrough) image; its µm/pixel is read from the 2 µm `scalefactors_json.json`. See [Image Modalities](../run_together_images.md).
- Square (8/16 µm) and segmented-cell imports are wired by the profile and appended to the catalog when the directories exist.

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
- Step-by-step (manual modules): [Visium HD end-to-end tutorial](../../vignettes/pipelines/visiumhd.md)
