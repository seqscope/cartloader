# Image Modalities in `run_together`

Morphology / histology images (DAPI, H&E, protein stains, …) are packaged as georeferenced PMTiles layers and added to each sample's `catalog.yaml`. This page is the **platform-agnostic** reference for how to declare them; the per-platform pages show the exact image files each platform ships.

Images are declared in one of three equivalent ways:

| Way | Flag | Use when |
|-----|------|----------|
| **Image sheet** | `--images images.tsv` | Multi-sample, or several images; one row per image. |
| **Inline CLI** | `--image key=value,…` (repeatable) | A single sample, a few images. |
| **Config JSON** | `"images": [ … ]` in `--config` | Alongside other JSON settings. |

All three share the **same field vocabulary**. Built-in profiles also **auto-detect** standard modalities from `in_dir` (e.g. Xenium's `morphology_focus_000{0..3}.ome.tif`), so you often need none of these — see the per-platform pages.

---
## Fields

| Field (aliases) | Required | Meaning |
|---|:---:|---|
| `type` | ✅ | Modality — sets the **default color** and **kind** from the registry (below). |
| `source` / `src` / `tif` | ✅ | Image path, `.tif`/`.ome.tif` or `.png`. Relative paths resolve against the sample's `in_dir`. |
| `sample` / `sample_id` | (sheet only) | Which sample the image belongs to; `*` or blank = every sample. |
| `id` / `img_id` | | Catalog key / output basename. **Defaults to `type`.** |
| `color` | | Hex color, **overrides** the type's default (colorized `single` kinds only). |
| `kind` | | `single` \| `rgb` \| `prebuilt`. **Defaults to the type's registry kind.** |
| `transform` (or the platform column, e.g. `merfish_csv`) | | Geometric transform file → the profile's import flag (e.g. MERSCOPE `--micron2pixel-csv`). |
| `shrink_factor` | | Downscale factor for very large mosaics (else the profile's `image_defaults`). |
| `high_memory` | | `true` to allow a high-memory import path for big images. |
| `convert` | | `ome2png` \| `png2pmtiles` \| `none`. Defaults from the source extension (`.png` → none; `.tif` → `ome2png`). |
| `um_per_pixel` | | Microns per pixel, for a plain (non-OME) image that carries no pixel-size metadata. Sets `--px-per-um-x/y` (single-channel) or `--um-per-pixel` (rgb). E.g. `0.5` for a Stereo-seq `*_regist.tif`. |
| `georef_plain` | | `true` to map the upper-left corner to (0,0) and the lower-right to the image size, scaled by `um_per_pixel` (`rgb` images). |

Blank cells (`` / `-` / `.` / `NA`) mean *unset*.

### Kinds

| `kind` | Behavior | Typical use |
|--------|----------|-------------|
| `single` | grayscale → colorize with `color` | DAPI, protein/RNA stains (OME-TIFF) |
| `rgb` | passthrough (no colorize), georeference | H&E histology |
| `prebuilt` | copy an existing `.pmtiles` straight into the catalog | pre-rendered layers |

### Type registry (default color + kind)

`type` looks up a default color and kind in `assets/run_together_image_types.json` — edit that file to add or adjust types (keys starting with `_` are ignored):

| `type` | default color | kind |
|--------|:---:|:---:|
| `dapi` | `0F73E6` | single |
| `boundary` | `F300A5` | single |
| `rna` | `A4A400` | single |
| `protein` | `008A00` | single |
| `cellbound1` / `cellbound2` / `cellbound3` | `009E73` / `E69F00` / `CC79A7` | single |
| `polyt` | `00FFFF` | single |
| `hne` | — (RGB) | rgb |

An **unregistered** `type` with no explicit `color` is an error (no silent wrong color) — either add a `color` field or register the type.

---
## Image sheet (`--images images.tsv`)

One row per image, attached to samples by the `sample` column. This keeps the sample sheet from blowing up while staying flexible.

```
sample    type      source            color     merfish_csv
s1        dapi      /p/s1_dapi.tif     -         /p/s1_transform.csv
s1        protein   /p/s1_cd3.png      008A00    /p/s1_transform.csv
s2        hne       /p/s2_he.tif       -         -
```

- `sample = *` (or blank) applies a row to **every** sample.
- The transform column is named by the platform profile (`merfish_csv` for MERSCOPE); a generic `transform` column also works everywhere.

---
## Inline CLI (`--image`)

For a single sample, give each image as comma-separated `key=value` pairs; repeat `--image` per image (same fields as a sheet row, no `sample` needed):

```bash
cartloader run_together --platform merfish \
  --in-transcript /p/molecules.csv --in-cell-boundary /p/polys.csv \
  --image type=dapi,source=/p/dapi.tif,merfish_csv=/p/transform.csv \
  --image type=protein,source=/p/protein.png,color=008A00,transform=/p/transform.csv \
  --id MYSAMPLE --out-dir OUT
```

---
## Config JSON (`"images"`)

The same records as a list under `--config`. `source` is an explicit path (or `in_dir`-relative); `match` is an `in_dir`-relative glob used for auto-detection (first match whose file exists wins — same-`id` fallback entries are allowed):

```jsonc
{ "id": "dapi", "match": "morphology_focus/..._0000.ome.tif", "kind": "single", "color": "0F73E6", "convert": "ome2png" }
{ "id": "hne",  "source": "he.tif", "kind": "rgb", "convert": "png2pmtiles" }
{ "id": "boundary", "source": "prebuilt/boundary.pmtiles", "kind": "prebuilt" }
```

Per-sample / JSON entries **append** to (or override, by `id`) the profile's auto-detected images.

---
## What a row runs

A colorized `single` image runs roughly:

```
cartloader import_image --ome2png --png2pmtiles --georeference \
    --colorize <color> [--<transform-flag> <transform>] [--shrink-factor N] [--high-memory] \
    --in-img <source> --out-dir cartl/<sample>/ --img-id <id>
```

`shrink_factor` / `high_memory` fall back to the profile's `image_defaults` (e.g. MERSCOPE sets `shrink_factor=5.0`, `high_memory=true` for its large mosaics). An `rgb` image goes through the `image_png2pmtiles` (geotiff → mbtiles → pmtiles) path instead; a `prebuilt` image is copied straight in. See [`import_image`](./import_image.md).

---
## See also

- [Specifying Inputs](./run_together_inputs.md) — the three input modes and the sample sheet.
- Per-platform image files: [Xenium](./platforms/xenium.md), [Visium HD](./platforms/visium_hd.md), [CosMx SMI](./platforms/cosmx_smi.md), [MERSCOPE](./platforms/merscope.md).
- [`import_image`](./import_image.md) — the underlying module.
