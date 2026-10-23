# Image Modalities in `run_together`

Morphology / histology images (DAPI, H&E, protein stains, …) are packaged as georeferenced PMTiles layers and added to each sample's `catalog.yaml` in the `images` stage. This page is the **platform-agnostic** reference for how to declare them; the per-platform pages show the image files each platform ships.

Images are declared in one of three ways:

| Way | Flag | Use when |
|-----|------|----------|
| **Image sheet** | `--images images.tsv` | Multi-sample, or several images; one row per image. |
| **Inline CLI** | `--image key=value,…` (repeatable) | A single sample, a few images. |
| **Config JSON** | `"images": [ … ]` in `--config` (top level, or inside a sample entry) | Alongside other JSON settings; also how profiles declare auto-detected images. |

The sheet and `--image` share the **same field vocabulary** (below). Config JSON entries are used **as written**, with a few differences — see [Config JSON](#config-json-images). Built-in profiles also **auto-detect** standard modalities (e.g. Xenium's `morphology_focus_000{0..3}.ome.tif`, Stereo-seq's `*_regist.tif`), so you often need none of these.

| Platform | Images by default |
|----------|-------------------|
| [Xenium](./platforms/xenium.md) | auto-detected `morphology_focus` channels, or `morphology_focus.ome.tif` / `morphology.ome.tif` |
| [Visium HD](./platforms/visium_hd.md) | H&E from the sample's `hne` column / key, scaled from `scalefactors_json.json` |
| [MERSCOPE](./platforms/merscope.md) | none; declare them with the `merfish_csv` transform column |
| [Stereo-seq](./platforms/stereoseq.md) | auto-detected `{prefix}_HE/ssDNA/DAPI_regist.tif` at 0.5 µm/pixel |
| [Seq-Scope](./platforms/seqscope.md) | none; an H&E given explicitly is treated as already georeferenced |
| [CosMx SMI](./platforms/cosmx_smi.md) · [Illumina StrataMap](./platforms/illumina.md) · [Generic](./platforms/generic.md) | none; attach explicitly |

---
## Fields

Fields of an `--images` row or an `--image` value:

| Field (aliases) | Required | Meaning |
|---|:---:|---|
| `type` | ✅ | Modality — sets the **default color** and **kind** from the registry (below). |
| `source` / `src` / `tif` | ✅ | Image path. Relative paths are resolved against the sample's `in_dir` when not found as given. |
| `sample` / `sample_id` | (sheet only) | Which sample the image belongs to; `*` or blank = every sample. |
| `id` / `img_id` | | Catalog key / output basename. **Defaults to `type`.** |
| `color` | | Hex color, **overrides** the type's default (`single` kind only). |
| `kind` | | `single` \| `rgb` \| `prebuilt`. **Defaults to the type's registry kind.** |
| `transform` (or the platform column, e.g. `merfish_csv`) | | Geometric transform file → the profile's import flag (MERSCOPE: `--micron2pixel-csv`). `single` kind only. |
| `convert` | | `single` kind: `ome2png` converts the TIFF and supplies the georeferencing bounds. Defaults to `ome2png` for a `.tif`, `none` for a `.png`. Any value other than `ome2png` skips that step (see [Georeferencing](#georeferencing)). |
| `um_per_pixel` | | Microns per pixel, for a plain (non-OME) TIFF that carries no pixel-size metadata, e.g. `0.5` for a Stereo-seq `*_regist.tif`. `single`: sets `--px-per-um-x/y`. `rgb`: sets `--um-per-pixel`, and takes effect only together with `georef_plain`. |
| `shrink_factor` | | `single`: downscale factor for very large mosaics (else the profile's `image_defaults`). |
| `high_memory` | | `single`: `true` for the high-memory import path. Defaults to `true` whenever `use_middle_page` is in effect; `false` opts out. |
| `use_middle_page` | | `single`: `true`/`false`, pick the middle page of a multi-page (z-stack) OME-TIFF. **Defaults to `true` for files whose name ends in `morphology.ome.tif`** (Xenium's 3D DAPI stack, also under a prefix such as `GSM123_rep1_morphology.ome.tif`; `import_image` otherwise rejects it for lacking a page), `false` for everything else. |
| `georeferenced` | | `rgb`: `true` if the TIFF already carries the intended geotransform (e.g. a Seq-Scope H&E registered upstream): its corner coordinates are kept and only the CRS is assigned (`--georef-detect gtiff`). Else the profile's `image_defaults`. |
| `georef_detect` | | `rgb`: read the bounds from the file — `ome` (OME `PhysicalSizeX/Y`) or `gtiff` (as `georeferenced`). Else the profile's `image_defaults`. An `.ome.tif` defaults to `ome`. |
| `georef_plain` | | `rgb`: map the upper-left corner to (0,0) and the lower-right to the image size × `um_per_pixel`. Any non-empty value turns it on. |
| `swapxy` / `rotate` / `flip_vertical` / `flip_horizontal` | | `rgb`: orientation (see [Orientation](#orientation-rgb-images)). Else the profile's `image_defaults`. |
| `rescale` | | `rgb`: rescale mode for **16-bit** imagery (`auto` \| `linear` \| `log` \| `none`), forwarded to `geotiff2pmtiles`, which rejects 16-bit input without a range. Use `linear` with the range below. |
| `rescale_range` / `rescale_min` + `rescale_max` | | `rgb`: input value range for `rescale` (e.g. `0,65535`). Give `rescale_range` in a `--images` TSV or config; on a `--image` CLI value use `rescale_min`/`rescale_max` instead, since the range's comma would split the value. |
| `nodata` | | `rgb`: pixels whose bands all equal this value are written transparent (forwarded to `geotiff2pmtiles --nodata`), e.g. `0` for a black background. `none` disables a type default. Else the profile's `image_defaults`. |

Blank cells (`` / `-` / `.` / `NA`) mean *unset*. Fields marked `single` or `rgb` are accepted for the other kind but have no effect there.

### Kinds

| `kind` | Behavior | Typical use |
|--------|----------|-------------|
| `single` | grayscale → colorize with `color` (`import_image`) | DAPI, protein/RNA stains (TIFF / OME-TIFF) |
| `rgb` | passthrough, no colorizing (`image_png2pmtiles`) | H&E histology |
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
| `ssdna` | `F0F0F0` | single |
| `hne` | — (RGB) | rgb |
| `sbcd` | — (RGB) | rgb |

A registry entry may also carry `defaults`: per-image fields filled in for that type (an explicit field on the row wins, and they take precedence over the profile's `image_defaults`). `sbcd` (the Seq-Scope `spatula draw-3way` spatial-barcode PNG) uses them to declare a plain, palette-indexed PNG at 1 µm/pixel from the origin with a black background: `georeferenced=false`, `georef_plain=true`, `um_per_pixel=1`, `nodata=0`. See [Seq-Scope](./platforms/seqscope.md#spatial-barcode-image-typesbcd).

An **unregistered** `type` with no explicit `color` is an error (no silent wrong color) — either add a `color` field or register the type. An unregistered type defaults to `kind` `single`.

---
## Georeferencing

Every layer needs its bounds in transcript µm. Where they come from depends on the kind:

**`single`** — the OME→PNG step (`convert: ome2png`) supplies them: from the OME-TIFF's pixel-size metadata, from `um_per_pixel` for a plain TIFF, and through the platform transform (`merfish_csv`) when given. A source that skips this step (a `.png`, or `convert` set to anything else) reaches `image_png2pmtiles --georeference` with no bounds and the import fails, so single-channel layers need a TIFF source (or, in JSON, `extra_flags` that supply a `--georef-*` option).

**`rgb`** — the first that applies:

1. `georeferenced` is true → `--georef-detect gtiff` (keep the file's geotransform, assign the CRS).
2. `georef_detect` is set → `--georef-detect <value>`.
3. The source is an `.ome.tif` / `.ome.tiff` and none of `georef_plain`, `um_per_pixel`, `um_per_pixel_json` is set → `--georef-detect ome`.
4. `georef_plain` → `--georef-plain`, scaled by `um_per_pixel_json` (a JSON under `in_dir`, key `um_per_pixel_key`, default `microns_per_pixel`; read with `jq` — the Visium HD H&E setting) or by `um_per_pixel` (per image; `image_defaults.um_per_pixel` applies to `single` images only).

With none of these, `image_png2pmtiles` stops with "georeferencing bounds not found or ambiguous".

### Orientation (`rgb` images)

`rotate` (`90`, `180`, `270`, clockwise), `flip_vertical`, and `flip_horizontal` are forwarded to `image_png2pmtiles` (rotation precedes flips). `swapxy=true` transposes the image — for histology stored with its axes swapped relative to the transcript coordinates — and is shorthand for `rotate=90` plus `flip_vertical`, so it cannot be combined with an explicit `rotate` / `flip_*` (an error). Each falls back to the profile's `image_defaults`. They have no effect on `single` images. See [Image orientation](../faq/background_image_orientation.md) for how to choose.

---
## Image sheet (`--images images.tsv`)

One row per image, attached to samples by the `sample` column. This keeps the sample sheet from blowing up while staying flexible.

```
sample    type      source              color     merfish_csv
s1        dapi      /p/s1_dapi.tif      -         /p/s1_transform.csv
s1        protein   /p/s1_cd3.tif       008A00    /p/s1_transform.csv
s2        hne       /p/s2_he.ome.tif    -         -
```

- `sample = *` (or blank) applies a row to **every** sample.
- A single DAPI or H&E per sample can instead go straight into the **sample sheet** as a `dapi` / `hne` column (path only); see [Specifying Inputs](./run_together_inputs.md#the-sample-sheet-is-a-wide-table-of-input-roles) for how those are imported.
- The transform column is named by the platform profile (`merfish_csv` for MERSCOPE); a generic `transform` column also works on such a platform. On a platform without an `image_transform`, the column is ignored.

---
## Inline CLI (`--image`)

For a single sample, give each image as comma-separated `key=value` pairs; repeat `--image` per image (same fields as a sheet row, no `sample` needed):

```bash
cartloader run_together --platform merfish \
  --raw-transcript /p/detected_transcripts.csv --in-cell-boundary /p/cell_boundaries.csv \
  --image type=dapi,source=/p/mosaic_DAPI_z3.tif,merfish_csv=/p/micron_to_mosaic_pixel_transform.csv \
  --image type=protein,source=/p/mosaic_CD3_z3.tif,color=008A00,transform=/p/micron_to_mosaic_pixel_transform.csv \
  --id MYSAMPLE --out-dir OUT
```

---
## Config JSON (`"images"`)

The same records as a list under `--config` — at the top level (every sample) or inside a sample entry (that sample only). Profiles declare their auto-detected images the same way.

```jsonc
{ "id": "dapi",     "match": "morphology_focus/morphology_focus_0000.ome.tif", "kind": "single", "color": "0F73E6" }
{ "id": "ssdna",    "suffix": "_ssDNA_regist.tif", "kind": "single", "color": "F0F0F0", "um_per_pixel": 0.5 }
{ "id": "hne",      "source": "he.ome.tif", "kind": "rgb" }
{ "id": "boundary", "source": "prebuilt/boundary.pmtiles", "kind": "prebuilt" }
```

Each entry locates its file with exactly one of:

| Key | Meaning |
|-----|---------|
| `source` | Explicit path (relative → resolved against `in_dir` when not found as given). |
| `match` | A literal path under `in_dir` (not a glob). |
| `suffix` | Appended to the sample's `in_prefix` (Stereo-seq). |

JSON entries are **not** run through the type registry, so they differ from sheet rows:

- `id` is required (it does not default to `type`), and `type` sets nothing: give `kind` (default `single`) and, for `single`, `color` (required).
- `convert` defaults to `ome2png` whatever the extension.
- The `transform` / platform column is not recognized; pass a transform as `extra_flags` (e.g. `["--micron2pixel-csv", "/p/transform.csv"]`) or use a sheet row.
- `extra_flags` (raw `import_image` flags, `single` kind only) and `match` / `suffix` exist only here.

The other fields (`um_per_pixel`, `shrink_factor`, `high_memory`, `use_middle_page`, georeferencing, orientation, rescale) work as in the table above.

**Merging.** A top-level `images` entry replaces every profile entry with the same `id` and new ids are appended; `{ "replace": [ … ] }` discards the profile's list. See [List assembly](./run_together_inputs.md#list-assembly-rule-append-by-default-keyed-by-id).

---
## Which entry wins

For each sample the candidates are considered in this order: its own JSON `images`, then the `--images` / `--image` rows that target it, then its sample-sheet `dapi` / `hne` columns, then the profile/config `images`. For each `id`, the **first entry whose file exists** is used and the rest are ignored — so an explicit image overrides a profile's auto-detected one with the same id, and a profile can list several same-`id` fallbacks (Xenium's `dapi`).

!!! warning "A missing file is skipped silently"
    An entry whose file does not exist produces no layer and no error. If a layer is missing from the catalog, check the path (relative paths are tried as given, then under `in_dir`).

---
## What a row runs

A colorized `single` image runs:

```
cartloader import_image --ome2png [--skip-image-errors] --png2pmtiles --georeference \
    --in-img <source> --out-dir cartl/<dir>/ --img-id <id> \
    --upper-thres-quantile 0.95 --level 0 --colorize <color> --transparent-below 5 \
    [--micron2pixel-csv <transform>] [--px-per-um-x/-y <1/um_per_pixel>] \
    [--shrink-factor N] [--high-memory] [--use-middle-page] [<extra_flags>]
```

A file ending in `morphology.ome.tif` (Xenium's multi-page DAPI z-stack, with or without a prefix) gets `--use-middle-page` and `--high-memory`, whichever way it was supplied, unless `use_middle_page` / `high_memory` say otherwise. `shrink_factor`, `high_memory`, and `um_per_pixel` fall back to the profile's `image_defaults` (e.g. MERSCOPE sets `shrink_factor=5.0`, `high_memory=true` for its large mosaics; Stereo-seq sets `um_per_pixel=0.5`).

An `rgb` image runs `cartloader image_png2pmtiles --geotif2mbtiles --mbtiles2pmtiles --georeference` with the bounds, orientation, and rescale options above. A `prebuilt` image is copied in as `<id>.pmtiles`. Each layer is then added to `catalog.yaml` under its `id`. See [`import_image`](./import_image.md).

With **`--skip-image-errors`**, a `single` image whose OME-TIFF cannot be read during the OME→PNG step is skipped with a warning instead of failing the run, and any image whose PMTiles was not produced is left out of the catalog. Failures of `rgb` or `prebuilt` images still stop the run. Re-run with the flag to regenerate the Makefile and resume past a corrupt image.

---
## See also

- [Specifying Inputs](./run_together_inputs.md) — the three input modes and the sample sheet.
- Per-platform image files: [Xenium](./platforms/xenium.md), [Visium HD](./platforms/visium_hd.md), [MERSCOPE](./platforms/merscope.md), [Stereo-seq](./platforms/stereoseq.md), [Seq-Scope](./platforms/seqscope.md).
- [`import_image`](./import_image.md) — the underlying module.
