# Unified Pipeline Orchestrator (`run_together`)

## Overview

`run_together` is the **recommended way to run CartLoader**. One command takes raw platform output to a deployable CartoScope directory — **ingest → FICTURE → cell analyses → asset packaging → image layers → (optional) annotation and upload** — for one sample or many, on any [supported platform](./supported_platforms.md). It writes the whole run as **one Makefile** (`run_together.mk`) and executes it, so a run is **resumable** and **parallel**, and a multi-sample run trains **one joint FICTURE model**.

Use it whenever your data comes from a platform with a built-in profile, or can be described with the [`generic`](./platforms/generic.md) profile. Run the [individual commands](../vignettes/step_by_step.md) instead only when you need control over a stage that `run_together` does not expose.

A minimal single-sample run:

```bash
cartloader run_together --platform 10x_xenium \
    --in-dir /data/xenium/outs --out-dir /work/xenium-lung
```

The profile auto-detects every input inside `--in-dir` and supplies the FICTURE defaults (here `width=12`, `n_factor=12,24,48`). The deployable output is `/work/xenium-lung/cartl/xenium-lung-rep1/` (see [Output naming](#output-naming)).

- New to CartLoader? Start with [Get Started](../vignettes/quickstart/intro.md).
- Worked examples: [Single Sample with run_together](../vignettes/pipelines/run_together.md) · [Multiple Samples with run_together](../vignettes/pipelines/run_together_multi.md).

### Read next

| To… | See |
|-----|-----|
| Specify inputs (single-sample CLI, sample sheet, or config JSON) | [**Specifying Inputs**](./run_together_inputs.md) |
| Attach DAPI / H&E / protein images | [**Image Modalities**](./run_together_images.md) |
| See exact inputs and examples for your platform | [**Platforms**](./supported_platforms.md): [Xenium](./platforms/xenium.md) · [Visium HD](./platforms/visium_hd.md) · [CosMx SMI](./platforms/cosmx_smi.md) · [MERSCOPE](./platforms/merscope.md) · [Illumina StrataMap](./platforms/illumina.md) · [Stereo-seq](./platforms/stereoseq.md) · [Seq-Scope](./platforms/seqscope.md) · [Generic](./platforms/generic.md) |

---
## How a run is configured

> **Convenience for the default case; full control in JSON.**

Three layers assemble one canonical configuration, each refining the previous one:

1. **Profile** (`--platform <name>`, optionally deep-merged with `--platform-json`) — per-platform defaults: auto-detected inputs, default analyses, image conventions.
2. **CLI** — the common knobs (`--in-dir`/`--samples`, `--width`/`--n-factor` or `--project-models`, `--out-dir`, …).
3. **Config JSON** (`--config`) — augments or fully specifies anything.

The layers compose: set the base run on the CLI and add only extra analyses or images in JSON. The fully assembled result is written to `run_together.resolved.json` for provenance.

---
## Requirements

- `make`, `python3`, and `pigz` (every stage is run with `--gzip pigz`).
- The tools bundled as repository submodules — `spatula`, `punkst` (FICTURE2), `pmpoint`, `tippecanoe`, `pmtiles`, `geotiff2pmtiles`. The CartLoader modules that `run_together` calls look for them at their built locations under `submodules/`, so build them once with `submodules/build.sh` (see [Installation → 4](../installation.md#4-initializing-and-building-submodules)). `run_together` has no per-tool path flags.
- GDAL command-line tools (`gdal_translate`, `gdaladdo`, `gdalinfo`) for image layers; `jq` for the Visium HD H&E scale.
- Only when used: the SAW binary (`--saw`, Stereo-seq GEF input), the `aws` CLI (`--s3-upload`), and annotation API access (`--anno` / `--anno-deep`).

The Docker image `hyunminkang/cartloader:latest` (dated tags such as `hyunminkang/cartloader:20261007-1432` for reproducibility) contains all of the above. Its entry point is `cartloader`:

```bash
docker run --rm -v /data:/data -v /work:/work hyunminkang/cartloader:latest \
    run_together --platform 10x_xenium --in-dir /data/xenium/outs --out-dir /work/xenium-lung
```

---
## Stages, selection & resume

| Stage | What runs | Granularity |
|-------|-----------|-------------|
| `ingest` | `sge_convert` (or `reformat_cosmx`, or SAW + `sge_convert` for Stereo-seq), plus any `.parquet` → `.csv.gz` conversion, `spatula tsv-add-cell-id`, and `convert_cellxgene` | per sample |
| `ficture` | `run_ficture2_multi`, one call per FICTURE analysis, run sequentially | per output directory |
| `cells` | [`run_ficture2_multi_cells`](./run_ficture2_multi_cells.md), one call per active cell analysis | per output directory |
| `cartload` | `run_cartload2` (single sample) or [`run_cartload2_multi`](./run_cartload2_multi.md) (joint run) | per sample / per joint run |
| `images` | image layers (`import_image` / `image_png2pmtiles`), Visium HD square and cell imports, per-sample Xenium Ranger imports on joint runs | per sample |
| `anno` | [`anno_cartload_folder`](./anno_cartload_folder.md) (opt-in) | per sample / once per joint run |
| `upload` | `aws s3 cp` (opt-in) | per sample |

- **Resume:** re-run the same command (or `make -f OUT/run_together.mk -j N`). Completed steps are skipped via flag files in `OUT/mk/`.
- **Changed settings are not detected:** the flag files do not record parameters. After changing one, delete the affected `mk/*.done` files (and those of later stages), or use `--restart`.
- `--only ingest,ficture` / `--skip images` — run a subset. Excluded upstream stages are assumed done (their prerequisites are pruned so `make` does not fail). When `ingest` is not selected, the transcript-source check is skipped too, since a resumed run reuses the converted transcripts.
- `--restart` — rebuild everything (`make -B`).
- `--dry-run` — write the Makefile and print the commands (`make -n`) without executing. Planning still writes the list files under `OUT/tsv/`.

Packaging bundles everything produced: every FICTURE model's pixel decode plus every cell analysis that ran.

---
## Output naming

| Run | Sample id (when not given) | Packaged directory |
|-----|----------------------------|--------------------|
| **Single sample**, `--out-dir OUT/name` | `rep1` | `OUT/name/cartl/name-<id>/` — or `cartl/<id>/` when the out-dir basename equals the id |
| **Joint run** (≥ 2 samples sharing `--out-dir OUT/name`) | `in_dir` / `in_prefix` basename | `OUT/name/cartl/name-<id>/` per sample, plus `cartl/multi-catalog.yaml` |
| **`--out-root ROOT`** (independent per-sample models) | `in_dir` / `in_prefix` basename | `ROOT/<id>/cartl/<id>/` |

- A lone `--out-dir` sample without an id — from `--in-dir`, the `--in-*` file flags, a one-row sample sheet, or one config sample — is always `rep1`; set `--id` (or an `id` column/key) to name it. The descriptive name lives in `--out-dir`.
- The packaged directory name is also the catalog id (`run_cartload2 --id`).
- A **joint run** packages every sample with one [`run_cartload2_multi`](./run_cartload2_multi.md) call. Each sample is a self-contained `cartl/<multi_id>-<sample_id>/` (`multi_id` = out-dir basename), and `cartl/multi-catalog.yaml` links every per-sample `catalog.yaml` and holds the shared factor files (`post`/`rgb`/`de`/`info`/`umap`) at the `cartl/` root, so `cartl/` deploys as one unit. `run_cartload2_multi` discovers samples and shared assets from the [`ficture.multi.params.json`](./run_ficture2_multi.md#multi-sample-manifest) manifest that `run_ficture2_multi` writes.

---
## Publishing

Publishing is **opt-in** and **CLI-only** (there is no config block for it). Two independent actions; enable either or both. When both run, `upload` waits on `anno` (which edits `catalog.yaml`).

### Annotation (`--anno`, `--anno-deep`)

- **`--anno`** — AI-annotate each packaged directory with [`anno_cartload_folder`](./anno_cartload_folder.md) (one `annotate_bulk_de_with_ai` call per factor). **Requires `--tissue` and `--organism`.** Writes `<factor>-alias-ai.tsv` files and records them under each factor's `alias_ai` key.
- **`--anno-deep`** — annotate with [`annotate_factors_with_llm`](./annotate_factors_with_llm.md) instead: all factors of a model in **one** deep-reasoning LLM call, with alternative interpretations and an interactive HTML report (`<factor>-alias-ai.html`). Implies `--anno`. **Requires `--organism`**, and `--tissue` unless every sample has a `tissue` (sample-sheet column or JSON sample key); a sample without one falls back to `--tissue`. On a joint run the per-sample tissues are written to `OUT/mk/anno.tissues.tsv` and passed as `--sample-sheet`.

| Option | `--anno` default | `--anno-deep` default |
|--------|------------------|-----------------------|
| `--anno-api-type` | `umgpt` | `claude` |
| `--anno-model` | `claude-opus-4-7` | `claude-opus-5-5` |
| `--anno-effort` | — | `high` (`low`, `medium`, `high`, `xhigh`, `max`) |
| `--anno-threads` | `10` | `10` |

A joint run annotates the shared factors **once** at the `cartl/` root (`--multi-sample`) and reuses them in every sample. Factors that already carry `alias` (e.g. from a cell analysis `alias` file) or `alias_ai` are skipped.

### S3 upload (`--s3-upload`)

Each packaged sample directory is copied to

```
<s3-prefix>/batch=<batch>/<collection>/<packaged directory name>/
```

For a joint run, `multi-catalog.yaml` and the shared factor files are also uploaded to the parent `<s3-prefix>/batch=<batch>/<collection>/`, so its relative pointers resolve.

| Option | Default |
|--------|---------|
| `--s3-prefix` | `s3://cartostore/data` |
| `--batch` | current `YYYY_MM` |
| `--collection` | single-sample `--out-dir` run: the out-dir's **parent** basename · joint `--out-dir` run: the out-dir basename · `--out-root` run: the out-root basename |
| `--aws-profile` | `default` |
| `--aws` | `aws` |
| `--s3-jobs` | `4` (parallel copies per upload) |

For example, a single-sample run with `--out-dir /work/projA/lung` uploads `cartl/lung-rep1/` to `s3://cartostore/data/batch=<YYYY_MM>/projA/lung-rep1/`.

```bash
cartloader run_together --platform 10x_xenium --samples samples.tsv --out-dir OUT \
    --width 12 --n-factor 24 \
    --anno --tissue "Kidney" --organism human \
    --s3-upload --collection my-collection
```

---
## Command-line parameters

Grouped as in `cartloader run_together --help`. A CLI flag wins over the same setting in `--config` or the profile, except where noted.

### Run options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--dry-run` | off | Write the Makefile and print the commands without executing. |
| `--restart` | off | Ignore existing outputs and rebuild (`make -B`). |
| `-j`, `--n-jobs` | `1` | Parallel `make` jobs; also each command's `--n-jobs` unless the config sets `resources.n_jobs`. |
| `--threads` | `4` | Each command's `--threads`, unless the config sets `resources.threads`. |
| `--makefn` | `run_together.mk` | Master Makefile name. |
| `--only` | all | Comma-separated stages to run exclusively (`ingest,ficture,cells,cartload,images,anno,upload`). |
| `--skip` | none | Comma-separated stages to skip. |
| `--skip-image-errors` | off | Tolerate an unreadable/corrupt OME-TIFF in a single-channel image import: it is skipped with a warning and left out of the catalog instead of failing the run. See [Image Modalities](./run_together_images.md#what-a-row-runs). |
| `--allow-unknown-config-keys` | off | Ignore, instead of rejecting, unrecognized keys in `--config` (top-level keys, sample keys, closed-schema blocks, and sample roles the platform cannot consume). See [Specifying Inputs → Validation](./run_together_inputs.md#validation-of-config-keys). |

!!! note "Resources: the config wins"
    `resources.n_jobs` / `resources.threads` in `--config` take precedence over `-j` / `--threads` for the values passed to each command. `make` itself always runs with the CLI `-j`.

### Input / output

| Parameter | Description |
|-----------|-------------|
| `--platform` | Profile name: `10x_xenium`, `10x_visium_hd`, `cosmx_smi`, `merfish`, `illumina`, `stereoseq`, `seqscope`, `generic`. May instead come from the config's `platform`. |
| `--in-dir` | Single sample: input directory; roles are auto-detected inside it. |
| `--in-prefix` | Single sample: input path prefix; roles are auto-detected by suffix (Stereo-seq). |
| `--in-transcript` | Single sample: an **already-ingested** transcript TSV (`X, Y, gene, count[, cell_id]`); skips ingest. |
| `--raw-transcript` | Single sample: a **raw** transcript CSV/TSV (or `.parquet`) to ingest through `sge_convert`. Mutually exclusive with `--in-transcript`. Not for `stereoseq`. |
| `--in-cell-xy` | Single sample: cell centroid / metadata file (`xy` role). |
| `--in-cell-boundary` | Single sample: cell boundary polygon file (`boundaries` role). |
| `--in-cellxgene` | Single sample: cell-by-gene matrix CSV, converted to a MEX that drives cell clustering. |
| `--geojson` | Single sample: GeoJSON polygon(s) in µm; only transcripts inside are kept ([polygon filtering](./run_together_inputs.md#polygon-filtering)). For a sheet, use its `geojson` column. |
| `--samples` | TSV sample sheet, one row per sample. |
| `--images` | TSV of image specs, one row per image (attached by its `sample` column). |
| `--image` | One image as comma-separated `key=value` pairs; repeatable. |
| `--out-dir` | Output directory (single sample, or a joint run with one shared model). |
| `--out-root` | Output root; each sample gets `<out-root>/<id>/` and an independent model. |
| `--id` | Sample id of the single sample built from the flags above (default `rep1` with `--out-dir`). Ignored for `--samples` / config samples. |
| `--config` | Config JSON that augments the profile and CLI. |
| `--platform-json` | JSON deep-merged over the built-in profile named by `--platform` (see [Supported Platforms](./supported_platforms.md#custom-profiles)). |
| `--saw` | Path to the SAW binary; required only when a Stereo-seq `.gef` is converted. |

`--in-dir`, `--in-prefix`, and the `--in-*`/`--raw-transcript` flags together describe **one** sample. It is added to any `--samples` rows and config `samples`; all of them sharing one `--out-dir` form a joint run. Details: [Specifying Inputs](./run_together_inputs.md).

### Ready-made cell input lists

Each supplies a `run_ficture2_multi_cells` list file verbatim instead of the one `run_together` derives from the samples, for every cell analysis. See [Specifying Inputs → Supplying your own `--list-*` files](./run_together_inputs.md#supplying-your-own-list-files).

| Parameter | Feeds | Line format |
|-----------|-------|-------------|
| `--list-cluster` | `--list-cluster` | `SAMPLE_ID<TAB>CLUSTER_FILE` (external cluster labels; without it clusters are computed on demand) |
| `--list-xy` | `--list-xy` | `SAMPLE_ID<TAB>XY_FILE` |
| `--list-boundaries` | `--list-boundaries` | `SAMPLE_ID<TAB>BOUNDARY_FILE` |
| `--list-mex` | `--mex-list` | `SAMPLE_ID<TAB>MEX_DIR` (or a bcd/ftr/mtx triple) |
| `--list-cell-tsv` | `--tsv-list` | `SAMPLE_ID<TAB>CELL_TSV` |

### Input column overrides

They set profile-level values, so they apply to **every** sample in the run.

| Parameter | Sets |
|-----------|------|
| `--colname-transcript-x` / `-y` / `-feature` / `-count` | Raw transcript columns read by `sge_convert` (`ingest.csv_colnames`). Without a count column each row counts as 1. |
| `--colname-transcript-cell` | A cell-id column already present in the raw transcript: carried to transcript column 5 and used for cell analysis without `tsv-add-cell-id` (`ingest.csv_colname_cell`). |
| `--colname-xy-cell` / `-x` / `-y` | Columns of the `xy` file (`roles.xy`); `--colname-xy-cell ''` for an unnamed index column. |
| `--colname-boundary-cell` / `-x` / `-y` | Columns of the `boundaries` file read by `spatula tsv-add-cell-id` (`roles.boundaries`). |

### FICTURE mode

| Parameter | Description |
|-----------|-------------|
| `--width` | De-novo hexagon width(s) in µm, comma-separated (overrides the profile's training entry). |
| `--n-factor` | De-novo factor count(s), comma-separated; widths × factors are all trained. |
| `--project-models` | Projection only: existing FICTURE dir(s), comma-separated; every model in each `ficture.params.json` is reused. No training. |
| `--no-ficture` | Tiling only: package transcripts, raster, and images with no factor layers; cell analyses are skipped. Not combinable with `--n-factor` / `--project-models`. |
| `--segment-10x` | Also export the hexagon files as 10x MEX directories (any mode). |
| `--segment-width-10x` | Widths to export with `--segment-10x` (default: each analysis' own width). |
| `--tile-size` | punkst tiling size in µm (default: `run_ficture2_multi`'s `500`). |
| `--tile-buffer` | punkst per-tile line buffer (default: `run_ficture2_multi`'s `1000`). |

See [Specifying Inputs → FICTURE mode](./run_together_inputs.md#ficture-mode-de-novo-vs-projection).

### Common decode overrides

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--exclude-feature-regex` | <code>^(Gm[0-9]&#124;MT-&#124;mt-&#124;Rps&#124;Rpl)</code> | Features kept out of the factor models (data and packaging keep them). `''` disables. |
| `--include-feature-list` / `--exclude-feature-list` | none | File of feature names (one per line) to restrict the factor models to / exclude from them. |
| `--min-ct-per-unit-hexagon` | `50` | Minimum count per hexagon (all genes). |
| `--min-ct-per-unit-train` | punkst's own | Minimum count per hexagon during training, over the features the model uses. |
| `--cell-min-cell-count` | `run_ficture2_multi_cells`' own (`50`) | Minimum count per cell for cells built from a transcript or cell TSV. |
| `--cell-min-feature-count` | `run_ficture2_multi_cells`' own (`20`) | Minimum total count per feature in the cell analysis. |
| `--always-single-molecule` / `--never-single-molecule` | ON for pixel FICTURE, OFF for cells | Force single-molecule mode on/off for both. |
| `--colname-cell` | none | `generic` platform: turns on cell analysis from the transcript's 5th column (`cell_id`); also the default cell-id column of the `xy` file. See [Generic](./platforms/generic.md). |

### Ingest feature filters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ingest-exclude-feature-regex` | <code>^(Unassigned&#124;Neg&#124;BLANK&#124;Blank&#124;Intergenic&#124;Deprecated&#124;System&#124;NCS-&#124;NCP-)</code> | Features dropped while the transcript TSV is written (gone from every later stage). Replaces `sge_convert`'s per-platform default. `''` disables. |
| `--ingest-include-feature-regex` | none | Features kept while the transcript TSV is written. |
| `--ingest-include-feature-list` / `--ingest-exclude-feature-list` | none | Files of feature names to keep / drop at ingest. |

Ingest filters (and the two coordinate options below) are ignored, with a warning, for a sample that supplies an already-ingested transcript. See [Specifying Inputs → Feature filtering](./run_together_inputs.md#feature-filtering-two-independent-layers).

### Ingest coordinate options

| Parameter | Description |
|-----------|-------------|
| `--jitter-xy` | Uniform random ± µm offset added to each transcript's X and Y at ingest (default `0`, off). Not supported on `cosmx_smi`. See [Coordinate jitter](./run_together_inputs.md#coordinate-jitter-at-ingest). |
| `--units-per-um` | Input coordinate units per µm (default: each platform's convention). Not supported on `cosmx_smi`; no effect on `seqscope` (the profile's own `--units-per-um 1000` takes precedence) or `10x_visium_hd` (derived from the scale factors JSON). See [Coordinate units](./run_together_inputs.md#coordinate-units-at-ingest). |

### Packaging overrides

| Parameter | Description |
|-----------|-------------|
| `--bin-count` | Number of gene bins for the point PMTiles layers (`cartload.bin_count`; `500` in every built-in profile). |

### Publish

`--anno`, `--anno-deep`, `--tissue`, `--organism`, `--anno-api-type`, `--anno-model`, `--anno-effort`, `--anno-threads`, `--s3-upload`, `--collection`, `--batch`, `--s3-prefix`, `--aws-profile`, `--aws`, `--s3-jobs` — see [Publishing](#publishing).

---
## Output

```
out_dir/
├── run_together.mk               # the generated pipeline
├── run_together.resolved.json    # fully assembled config (provenance)
├── mk/                           # per-stage flag files (sge.<id>.done, ficture.done, cells.done, …)
├── tsv/                          # tsv/<id>/: per-sample ingest output; list files
│                                 #   (in_list.tsv, in_<role>.<analysis>.tsv, alias.<analysis>.tsv)
├── fic/                          # FICTURE results; per sample under fic/samples/<id>/
│   └── ficture.multi.params.json     # shared multi-sample manifest
└── cartl/                        # ← deploy this directory
    ├── <out_dir basename>-<id>/      # single sample: catalog.yaml + PMTiles
    ├── <multi_id>-<sample_id>/       # joint run: one self-contained dir per sample
    └── multi-catalog.yaml            # joint run only: links per-sample catalogs + shared assets
```

For `--out-root` runs this layout is created under `<out_root>/<id>/` for each sample, and the master Makefile and `run_together.resolved.json` are written to `<out_root>/`.

---
## See also

- [Specifying Inputs](./run_together_inputs.md) · [Image Modalities](./run_together_images.md) · [Supported Platforms](./supported_platforms.md)
- Guides: [single sample](../vignettes/pipelines/run_together.md) · [multiple samples](../vignettes/pipelines/run_together_multi.md)
- [Legacy orchestrators](./legacy.md) (`run_xenium`, `run_visiumhd`) and the [individual commands](../vignettes/step_by_step.md)
