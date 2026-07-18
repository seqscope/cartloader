# Supported Platforms & Inputs

This page is the canonical reference for **what `run_together` expects on disk** for each platform. The convenience of the one-command run comes from **platform profiles** that read a platform's standard output layout automatically — this table tells you what that layout must be.

For the profile mechanics (the canonical config, layer merging, custom profiles), see the [`run_together` reference](./run_together.md#the-canonical-configuration).

---
## Profile status

| `--platform` | `run_together` profile | Ingest via `sge_convert` |
|--------------|:----------------------:|:------------------------:|
| `10x_xenium` | ✅ built-in | ✅ |
| `10x_visium_hd` | ✅ built-in | ✅ |
| `cosmx_smi` | ✅ built-in | via `reformat_cosmx` |
| `merfish` (Vizgen MERSCOPE) | ✅ built-in | ✅ (`vizgen_merscope` preset) |
| `bgi_stereoseq` | ⏳ planned | ✅ |
| `seqscope` | ⏳ planned | ✅ |
| `pixel_seq`, `nova_st` | ⏳ planned | ✅ |
| `generic` (custom CSV/TSV) | via custom profile | ✅ |

- **built-in** — `run_together --platform <name>` works with no extra configuration.
- **planned** — `sge_convert` already supports ingest for these platforms; a `run_together` profile is not yet shipped. Use them today either through [`sge_convert`](./sge_convert.md) + the individual modules, or by supplying your own `--platform-json` file (see the [canonical config](./run_together.md#the-canonical-configuration)).

---
## Built-in profiles

### `10x_xenium`

Point `--in-dir` at the Xenium Ranger output directory. First match wins where several paths are listed; optional files are skipped when absent.

| Purpose | Path under `--in-dir` |
|---------|-----------------------|
| Transcripts | `transcripts.csv.gz`, or `transcripts.parquet`, or `transcripts/transcripts.parquet` |
| Cell boundaries | `cell_boundaries.csv.gz` |
| Cell centroids | `cells.csv.gz` (columns `x_centroid`, `y_centroid`) |
| Cluster labels | `analysis/clustering/gene_expression_graphclust/clusters.csv` → `xeniumranger` factor. Single-sample: decoded via `run_ficture2_multi_cells`. Joint run: these clusters are sample-specific, so each sample is imported per-sample via `import_xenium_cell` (sheet columns `xy`/`boundaries`/`clusters` are forwarded as `--csv-*` overrides, so GEO-style Ranger outputs with non-standard filenames also work). |
| Morphology images | `morphology_focus/morphology_focus_000{0,1,2,3}.ome.tif` → `dapi`/`boundary`/`rna`/`protein`; or a single `morphology_focus.ome.tif` / `morphology.ome.tif` → `dapi` |

Defaults: FICTURE `width=12`, `n_factor=12,24,48`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`.

### `10x_visium_hd`

Point `--in-dir` at the Space Ranger `outs/` directory. Coordinates are scaled by 2 automatically (`--scale-xy 2.0` at ingest, `--sge-scale 2` at packaging).

| Purpose | Path under `--in-dir` |
|---------|-----------------------|
| Scale factors | `binned_outputs/square_002um/spatial/scalefactors_json.json` |
| Count matrix | `binned_outputs/square_002um/filtered_feature_bc_matrix/` |
| Bin positions | `binned_outputs/square_002um/spatial/tissue_positions.parquet` |
| Square layers | `binned_outputs/square_008um`, `binned_outputs/square_016um` (imported if present) |
| Segmented cells | `segmented_outputs/` (imported if present) |
| H&E image | Provided per sample via the config `hne` field; µm/pixel read from the 2 µm `scalefactors_json.json` |

Defaults: FICTURE `width=12`, `n_factor=24,48,96`, `decode_scale=2`; packaging with `--use-pmpoint --bin-count 500 --sge-scale 2`.

### `cosmx_smi`

Point `--in-dir` at the AtoMx / CosMx SMI flat-file export directory. Unlike the other platforms, ingest does not run `sge_convert`: a dedicated [`reformat_cosmx`](./reformat_cosmx.md) step reads the three raw CSVs and, converting global-pixel coordinates to microns (`0.12028 µm/px`), writes the transcript TSV plus the cell-metadata (xy) and polygon (boundaries) files that FICTURE and cell decode consume. Filenames are matched by glob (first match wins).

| Purpose | Default glob under `--in-dir` | reformat flag |
|---------|-----------------------|:---:|
| Transcripts | `*tx_file.csv.gz` (columns `fov`, `cell_ID`, `target`, `x_global_px`, `y_global_px`, `z`) → `*.transcripts.tsv.gz` | `--tx` |
| Cell centroids | `*metadata*.csv.gz` (columns `fov`, `cell_ID`, `CenterX_global_px`, `CenterY_global_px`) → `*.metadata.csv.gz` (xy role, columns `X`/`Y`) | `--meta` |
| Cell boundaries | `*polygons.csv.gz` (columns `fov`, `cellID`, `x_global_px`, `y_global_px`) → `*.polygons.csv.gz` | `--poly` |

Cell ids are formed as `<fov>_<cell_ID>`; transcripts with `cell_ID = 0` are kept but tagged `UNASSIGNED`, and `System*` control probes are dropped at ingest. The `cartloader` cell factor is decoded from the metadata + polygon files via `run_ficture2_multi_cells`. For a joint multi-sample run of already-reformatted samples, supply the `transcript`/`cell_xy`/`cell_boundary` columns in the sample sheet (pointing at the `*.transcripts.tsv.gz` / `*.metadata.csv.gz` / `*.polygons.csv.gz` files) to skip re-ingest.

#### Overriding the input file patterns

CosMx exports don't always use the standard suffixes. Each input is matched by a **glob pattern** (exactly one file must match — zero or multiple is an error that lists what's present), and any pattern can be overridden with a small JSON config whose `ingest.inputs` block **deep-merges** over the profile — you only restate the patterns that differ:

```json
{
  "platform": "cosmx_smi",
  "ingest": { "inputs": { "--tx": "*tx_unique.csv.gz" } }
}
```

```bash
cartloader run_together --config cosmx_touchstone.json --in-dir IN --out-dir OUT
```

For example, an export with `..._tx.csv.gz`, `..._tx_unique.csv.gz`, `..._metadata.csv.gz`, and `..._polygons.csv.gz` needs only the `--tx` override above to pick `tx_unique` (the default `--meta`/`--poly` globs already match `_metadata`/`_polygons`); override `--meta`/`--poly` the same way when their names differ.

Defaults: FICTURE `width=12`, `n_factor=12,24,48`, `min_ct_per_unit_hexagon=100`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`.

### `merfish`

For MERFISH / Vizgen MERSCOPE. The inputs are **individual files that may be arbitrarily named and need not share a directory**, so point at each one explicitly (a single sample, no sample sheet needed):

```bash
cartloader run_together --platform merfish \
  --in-transcript   path/to/transcripts.csv \
  --in-cell-xy      path/to/cell_metadata.csv \
  --in-cell-boundary path/to/cell_boundaries.csv \
  --id  MYSAMPLE --out-dir OUT
```

| Flag | File | Expected columns (defaults) |
|------|------|-----------------------------|
| `--in-transcript` (required) | Raw transcript CSV, ingested through `sge_convert` (the `vizgen_merscope` preset) | `global_x`, `global_y`, `gene` (`global_z` optional) |
| `--in-cell-xy` (optional) | Cell metadata / centroids → xy role | `center_x`, `center_y`; the unnamed pandas-index first column holds the cell id |
| `--in-cell-boundary` (optional) | Cell boundary polygons → boundaries role | `cell_id`, `vertex_x`, `vertex_y` |

**Column names may be overridden per file** when they differ from the defaults above: `--colname-transcript-x/-y/-feature/-count`, `--colname-xy-cell/-x/-y` (use `--colname-xy-cell ''` for an unnamed index column), and `--colname-boundary-cell/-x/-y`.

Only `--in-transcript` is required; with just it, the run is pixel-level (FICTURE → packaging). Supplying `--in-cell-boundary` enables cell-level decode: unlike CosMx, MERFISH transcripts carry **no per-transcript cell assignment**, so ingest runs [`spatula tsv-add-cell-id`](../../submodules/spatula) after `sge_convert` to assign each transcript to the polygon that contains it (point-in-polygon; transcripts inside no cell are tagged `UNASSIGNED`), appending a `cell_id` column so the standard `cartloader` cell factor decodes via `run_ficture2_multi_cells` — exactly as for CosMx. `--in-cell-xy` adds the per-cell spatial cluster scatter. Coordinates are already in microns and are not shifted, so the three files stay aligned. Set nearest-cell fallback assignment (assign an outside-all-cells transcript to the closest cell within a distance) with `{ "ingest": { "assign_cell_id": { "expand_um": <µm> } } }`.

A standard MERSCOPE export directory also works via `--in-dir OUT`, which auto-detects `detected_transcripts.csv[.gz]`, `cell_metadata.csv`, and `cell_boundaries.csv`. For a joint multi-sample run, give `--samples` a sheet with an `in_dir` column per sample (each directory auto-detected as above); a per-row `cell_xy` / `cell_boundary` column overrides the auto-detected cell files. (Raw transcripts are auto-detected from `in_dir`, not passed as a sheet column — the sheet's `transcript` column means an *already-ingested* TSV and would skip `sge_convert`.)

Defaults: FICTURE `width=12`, `n_factor=12,24,48`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`.

#### Images

Morphology images also have non-standard paths, so they are given as a **separate image TSV** (`--images images.tsv`), one row per image, attached to samples by the `sample` column. This keeps the sample sheet from blowing up while staying flexible; the same records can equivalently be written as the config-JSON `images` list.

```
sample    type      source            merfish_csv           color
MYSAMPLE  dapi      /p/dapi.tif       /p/transform.csv      -
MYSAMPLE  protein   /p/protein.png    /p/transform.csv      008A00
```

| Column | Meaning |
|--------|---------|
| `sample` (required) | Sample id this image belongs to; `*` or blank = every sample |
| `type` (required) | Modality — sets the **default color** and kind. Registry: `dapi`→`0F73E6`, `boundary`→`F300A5`, `rna`→`A4A400`, `protein`→`008A00` (colorized single-channel); `hne` → multi-channel RGB passthrough (no color) |
| `source` / `src` (required) | Image path, `.tif` or `.png` (`--ome2png` is applied automatically for `.tif`, skipped for `.png`) |
| `id` | Catalog key / output basename; **defaults to `type`** |
| `color` | Hex, **overrides** the type's default. An unregistered `type` with no `color` is an error (no silent wrong color) |
| `merfish_csv` | The MERSCOPE `micron_to_mosaic_pixel_transform.csv` → `import_image --micron2pixel-csv`. The transform column name and flag are declared per platform by the profile's `image_transform` |

Each row runs `cartloader import_image --ome2png --png2pmtiles --georeference --colorize <color> --micron2pixel-csv <merfish_csv> …`. `--shrink-factor 5.0` and `--high-memory` are applied from the profile's `image_defaults` (large MERSCOPE mosaics); override per image with `shrink_factor` / `high_memory` columns.

For a single sample, the images can also be given **directly on the command line** with a repeatable `--image` — one image per flag, as comma-separated `key=value` pairs (the same fields as a TSV row; no `sample` needed). The transform accepts either the platform column name (`merfish_csv=`) or a generic `transform=`:

```bash
cartloader run_together --platform merfish \
  --in-transcript /p/molecules.csv --in-cell-boundary /p/polys.csv \
  --image type=dapi,source=/p/dapi.tif,merfish_csv=/p/transform.csv \
  --image type=protein,source=/p/protein.png,color=008A00,transform=/p/transform.csv \
  --id MYSAMPLE --out-dir OUT
```

---
## Platforms without a built-in profile yet

`sge_convert` converts all of the platforms in the status table to the unified transcript TSV that FICTURE and packaging consume — the per-platform column names, delimiters, and scaling are already encoded there (see the [`sge_convert` reference](./sge_convert.md)). Two ways to run them end-to-end today:

1. **Modules directly** — `sge_convert` → [`run_ficture2`](./run_ficture2.md) → [`run_cartload2`](./run_cartload2.md). See the platform starter tutorials for worked examples.
2. **A custom `run_together` profile** — copy a built-in profile (`assets/run_together_profiles/10x_xenium.json`), set `ingest.sge_platform` to the target platform and adjust the input paths / exclusion regex, then pass it with `--platform-json`. See the [canonical config](./run_together.md#the-canonical-configuration).

---
## See also

- [`run_together` reference](./run_together.md) — full options, tiers, and configuration schema.
- [Single-sample tutorial](../vignettes/pipelines/run_together.md) and [multi-sample tutorial](../vignettes/pipelines/run_together_multi.md).
- [`sge_convert` reference](./sge_convert.md) — per-platform ingest details.
