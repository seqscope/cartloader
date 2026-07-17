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
| `vizgen_merscope` | ⏳ planned | ✅ |
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

| Purpose | Path under `--in-dir` |
|---------|-----------------------|
| Transcripts | `*tx_file.csv.gz` (columns `fov`, `cell_ID`, `target`, `x_global_px`, `y_global_px`, `z`) → `*.transcripts.tsv.gz` |
| Cell centroids | `*metadata_file.csv.gz` (columns `fov`, `cell_ID`, `CenterX_global_px`, `CenterY_global_px`) → `*.metadata.csv.gz` (xy role, columns `X`/`Y`) |
| Cell boundaries | `*polygons.csv.gz` (columns `fov`, `cellID`, `x_global_px`, `y_global_px`) → `*.polygons.csv.gz` |

Cell ids are formed as `<fov>_<cell_ID>`; transcripts with `cell_ID = 0` are kept but tagged `UNASSIGNED`, and `System*` control probes are dropped at ingest. The `cartloader` cell factor is decoded from the metadata + polygon files via `run_ficture2_multi_cells`. For a joint multi-sample run of already-reformatted samples, supply the `transcript`/`cell_xy`/`cell_boundary` columns in the sample sheet (pointing at the `*.transcripts.tsv.gz` / `*.metadata.csv.gz` / `*.polygons.csv.gz` files) to skip re-ingest.

Defaults: FICTURE `width=12`, `n_factor=12,24,48`, `min_ct_per_unit_hexagon=100`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`.

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
