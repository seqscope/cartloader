# Multi-Sample Spatial Factor Inference Analysis with FICTURE

## Overview

`run_ficture2_multi` runs a joint FICTURE analysis across samples. It tiles each sample, builds joint hexagons (one or more widths), trains LDA models across samples per (width, factor count) — or projects a pre-trained model — then projects and pixel-decodes each sample, and writes per‑sample and shared JSON manifests.

Use it to learn shared spatial factors across samples and generate per‑sample outputs in a single, parallelizable Makefile. It also works with a single sample.

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) calls `run_ficture2_multi` for every FICTURE analysis (single- or multi-sample, de-novo or projection, and `--prepare-only` for `--no-ficture`). Use this page to run it directly in step-by-step or custom workflows.

!!! question
    - [What is `FICTURE` and `punkst`?](../faq/ficture.md)
    - [Why use Multi‑Sample FICTURE Analysis?](../faq/why_multi-sample.md)

---
## Requirements

- `--in-list`: a headerless, tab-delimited file with one line per sample: `<sample_id>` and the path to its transcript‑indexed SGE TSV (no blank or comment lines). Each transcript TSV has one header line; by default its first four columns are X, Y, feature, count (the [unified SGE format](./sge_convert.md#unified-sge-matrix)).
- `punkst`: `--ficture2` must point to a repository containing `bin/punkst` and `ext/py/factor_report.py` (default: the submodule `<cartloader_dir>/submodules/punkst`).
- `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`), `gzip`, `sort`, and `python3`.
- R (`Rscript` on the `PATH`) for the UMAPs, unless `--skip-umap` or `--prepare-only`.
- `cartloader` on the `PATH` (the Makefile calls `cartloader feature_select` and `cartloader write_json_for_ficture2_multi`).

---
## Actions

All actions run by default (UMAPs can be skipped with `--skip-umap`):

- **Tiling**: tiles each sample (`punkst multisample-prepare`) and writes the per-sample and shared feature lists.
- **Hexagons**: builds per-sample and joint hexagons at each `--width`.
- **Feature selection** (only with a [feature filter](#feature-filters)): resolves the filters into `multi.selected_features.tsv`.
- **LDA training**: trains a joint model for each (`--width`, `--n-factor`) pair, or projects `--pretrained-model` onto the joint hexagons; writes DE and a factor report.
- **Projection**: projects each model onto every sample's hexagons.
- **UMAP**: shared and per-sample UMAPs of the hexagon factor posteriors.
- **Decode**: pixel-level decoding of each sample with each model, with pseudobulk, DE, a factor report, and an image.
- **Manifests**: writes `samples/<sample>/ficture.params.json` and `ficture.multi.params.json`.
- **Hexagon MEX export** (optional, `--segment-10x`): converts each per-sample hexagon file into a 10x MEX directory.

`--prepare-only` stops after the hexagons: it writes each sample's `ficture.params.json` with only `in_sge` (tiled prefix, feature list, coordinate range) and an empty `train_params`. No model is trained, projected, or decoded, no UMAP is built, and feature filters are ignored; `--n-factor` is not required. This is what lets a dataset be **packaged without any factor analysis** — `run_cartload2` reads the tiled TSV as its molecule source, so it produces transcripts + raster + images and a catalog with no factor layers. Re-running without `--prepare-only` in the same `--out-dir` adds models later and re-uses the existing tiles/hexagons. See [`run_together --no-ficture`](./run_together_inputs.md#ficture-mode-de-novo-vs-projection).

`--segment-10x` additionally exports every per-sample hexagon file as a **10x MEX** directory, `samples/<sample>/<sample>.hex_<width>.mex/` (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`), via `spatula sptsv2mex`. The MEX is built from the very same `samples/<sample>/<sample>.hex_<width>.txt`/`.json` that the factor analysis uses — the hexagons are neither re-generated nor re-filtered, so the MEX contains exactly the hexagons that survive `--min-ct-per-unit-hexagon`. Barcodes are the hexagon centers as `x:y` (µm); features carry the gene name as both id and name, typed `Gene Expression`. By default every `--width` is exported; `--segment-width-10x` picks the widths instead (a width not in `--width` is simply added to the hexagon step, so its hexagon files get built too). The export depends only on the hexagon step, so it also works with `--prepare-only` — the way to get hexagon MEX files without any factor analysis. Each exported directory is recorded under `mex` (keyed by width) in the per-sample `ficture.params.json` and in `ficture.multi.params.json` (keyed by sample, then width).

### Reruns and done flags

Re-invoking the command with unchanged inputs leaves the tiles, hexagons, and every model downstream of them alone:

- Tiling is tracked by `multi.tiles.done`, which depends on `multi.tiles.params.json` (sample IDs and transcript paths, column indices, tile size/buffer) and on the transcript files, so a re-ingested sample is re-tiled.
- Each hexagon width is tracked by `multi.hex_<width>.done`, which depends on `multi.hex_<width>.params.json` (ordered sample list, width, `--min-count-per-sample`, `--min-ct-per-unit-hexagon`).
- Feature selection depends on `multi.selected_features.params.json` (the `feature_select` command) and on the list files.

A settings file is rewritten only when its content changes, so a changed setting, a changed sample list, or a new width reruns that step and everything built from it; a step whose outputs are missing reruns too. A changed sample list rebuilds every sample's tiles and hexagons. An output directory prepared by an older version (a single `multi.done` flag and no settings files) is adopted as up to date when its outputs exist and are not older than their inputs.

---
## Feature filters

`--include-feature-list`, `--exclude-feature-list`, `--include-feature-regex`, and `--exclude-feature-regex` restrict the features that **LDA training and projection** use:

- They are resolved by [`cartloader feature_select`](./feature_customization.md#resolving-filters-into-one-list-feature_select) (`--mode include`) against `multi.union_features.tsv` (every feature seen in any sample) into `multi.selected_features.tsv`, and the command is recorded in `multi.selected_features.params.json`.
- That list is passed as `lda4hex --features` to the joint training and to every per-sample projection; the pixel decode inherits the restriction through the model. Each sample's `ficture.params.json` records it as the FICTURE feature file.
- Tiling, the per-sample and shared feature lists, and the packaged genes keep **every** feature.
- Regexes use **substring** semantics (Python `re.search`); a list and a regex combine. Listed names absent from the data are ignored, and an empty selection is an error. Filters are ignored with `--prepare-only`.
- `--min-ct-per-unit-hexagon` counts all features, so when restricting to a small panel consider lowering `--min-ct-per-unit-train`, which counts only the selected features.

---
## Example Usage

De-novo training:

```bash
cartloader run_ficture2_multi \
  --in-list /path/to/samples.tsv \
  --out-dir /path/to/out/ficture_multi \
  --width 12 \
  --n-factor 12,24 \
  --exclude-feature-regex '^(Blank-|mt-|Gm[0-9])' \
  --threads 8 \
  --n-jobs 4
```

Projection of an existing model (no training):

```bash
cartloader run_ficture2_multi \
  --in-list /path/to/samples.tsv \
  --out-dir /path/to/out/ficture_multi \
  --width 12 \
  --pretrained-model /path/to/ref.model.tsv \
  --model-id ref
```

---
## Parameters

### Input/Output

- `--in-list` (str, required): Sample list (see [Requirements](#requirements)).
- `--out-dir` (str, required): Output directory.
- `--out-json` (str): Accepted but currently unused; the manifests are always `samples/<sample>/ficture.params.json` and `ficture.multi.params.json`.

### Tiling / Hexagons

- `--width` (str, required): Comma‑separated hexagon flat-to-flat widths in µm (e.g., `12,18`).
- `--colidx-x`, `--colidx-y`, `--colidx-feature`, `--colidx-count` (int, 1‑based): X, Y, feature and count columns of the transcript TSV (defaults: 1, 2, 3, 4).
- `--tile-size` (int): Tile size in µm. Must be at least 20× the side length (width / √3) of the widest hexagon, including `--segment-width-10x` (default: 500, raised automatically for wide hexagons to that minimum rounded up to 100, e.g. 600 for width 48).
- `--tile-buffer` (int): Per-tile buffer size (default: 1000).
- `--min-count-per-sample` (int): Minimum total count for a feature to enter the shared feature list `multi.features.tsv` (`punkst --min-total-count-per-sample`; default: 50).
- `--min-ct-per-unit-hexagon` (int): Minimum count per hexagon, over all features (default: 50).
- `--prepare-only` (flag): Tile and build hexagons only; write manifests with no models (`--n-factor` not required, UMAP forced off).
- `--segment-10x` (flag): Also export each per-sample hexagon file as a 10x MEX directory (`samples/<sample>/<sample>.hex_<width>.mex/`). Works with `--prepare-only`.
- `--segment-width-10x` (str): Comma‑separated widths (µm) to export with `--segment-10x` (default: the `--width` list); only valid with `--segment-10x`. Extra widths are added to the hexagon step.

### Training / Projection

- `--n-factor` (str): Comma‑separated factor counts (e.g., `12,24`). Required unless `--pretrained-model` or `--prepare-only`; not allowed with `--pretrained-model`. Every (width, factor count) pair is trained as model `t{width}_f{n_factor}`.
- `--pretrained-model` (str): A model TSV (feature by factor, header `Feature 0 1 ...`) to project instead of training; its factor count is read from the header.
- `--retrain` (flag): With `--pretrained-model`, use it as the prior and continue training instead of projecting only.
- `--model-id` (str): Name of the model in place of `t{width}_f{n_factor}`; only valid for a single model.
- `--min-ct-per-unit-train` (int): Minimum count per hexagon during training/projection, counted over the selected features (default: `punkst`'s default).
- `--minibatch-size` (int): Minibatch size (default: 500).
- `--train-epoch` (int): Number of training epochs (default: 2).
- `--seed` (int): Random seed (default: 1).
- `--cmap-file` (str): Color map TSV; its first `n_factor` colors are used (default: `<cartloader_dir>/assets/default_color_map.tsv`).
- `--skip-umap` (flag): Skip the UMAPs.
- Feature filters: `--include-feature-list`, `--exclude-feature-list`, `--include-feature-regex`, `--exclude-feature-regex` (see [Feature filters](#feature-filters)).

### Decoding / DE

- `--anchor-res` (int): Anchor resolution; the decode uses `width / anchor_res` sliding moves (default: 6).
- `--decode-pixel-res` (float): Pixel resolution (µm) of the decode (default: 0.5).
- `--single-molecule` (flag): Decode each molecule instead of pixels at `--decode-pixel-res`.
- `--decode-scale` (int): Scale from coordinates to pixels of the decode PNG (default: 1).
- `--de-min-ct-per-feature` (int): Minimum count per feature for differential expression (default: 20).
- `--de-max-pval` (float): P-value cutoff for differential expression (default: 1e-3).
- `--de-min-fold` (float): Fold-change cutoff for differential expression (default: 1.5).

### Run Options

- `--dry-run` (flag): Write the Makefile and print the commands only.
- `--restart` (flag): Ignore existing outputs and rerun all steps.
- `--threads` (int): Max threads per job (default: 8).
- `--n-jobs` (int): Parallel Makefile jobs (default: 2).
- `--makefn` (str): Makefile name under `--out-dir` (default: `run_ficture2_multi.mk`).

### Environment / Tools

- `--ficture2` (str): Path to the `punkst` repository (default: `<cartloader_dir>/submodules/punkst`).
- `--spatula` (str): Path to `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`).
- `--gzip` (str): Path to `gzip`; consider `pigz` for speed (default: `gzip`).
- `--sort` (str): Path to `sort` (default: `sort`).
- `--python` (str): Python 3 used for `punkst`'s factor report (default: `python3`).
- `--R` (str): R binary checked unless UMAPs are skipped (default: `R`); the UMAP scripts run with `Rscript` from the `PATH`.

??? note "Temporary files"
    `multisample-prepare` uses `<out-dir>/tmp`; `lda4hex` uses `<out-dir>/<model_id>_lda4hex` (joint) and `samples/<sample>/<sample>.<model_id>_lda4hex` (projection); `pixel-decode` uses `<out-dir>/tmp/<sample>_<model_id>_pixel`. The UMAP R scripts run with `TMPDIR` set to `<prefix>_create_umap`, `<prefix>_draw_umap` and `<prefix>_draw_umap_single`, where `<prefix>` is `<out-dir>/<model_id>` (shared) or `samples/<sample>/<sample>.<model_id>` (per sample). These directories are removed when their step finishes.

---
## Output

Outputs are written under `--out-dir`:

- Tiling and hexagons (joint):
    - `multi.union_features.tsv`: Every feature seen in any sample, with its total and per-sample counts.
    - `multi.features.tsv`: Shared feature list.
    - `multi.hex_<width>.txt` and `multi.hex_<width>.json`: Joint hexagons and metadata per width.
    - `multi.tiles.done`, `multi.tiles.params.json`, `multi.hex_<width>.done`, `multi.hex_<width>.params.json`: Done flags and their settings files (see [Reruns and done flags](#reruns-and-done-flags)).
    - `multi.selected_features.tsv` and `multi.selected_features.params.json`: With feature filters only.

- Per sample, in `samples/<sample>/`:
    - `<sample>.tiled.tsv`, `.tiled.index`, `.tiled.coord_range.tsv`, `.tiled.features.tsv`: Tiled transcripts (read by packaging).
    - `<sample>.hex_<width>.txt` and `.json`: The sample's hexagons.
    - `<sample>.hex_<width>.mex/`: 10x MEX export (`--segment-10x`).

- Joint LDA (per width × factor count), prefix `<model_id>` (default `t{width}_f{n_factor}`):
    - `.model.tsv`: Feature-by-factor weights.
    - `.results.tsv.gz`: Factor posteriors of the joint hexagons, with top factor columns.
    - `.cmap.tsv`: Factor colors.
    - `.bulk_chisq.tsv`: Per‑factor differential expression.
    - `.factor.info.tsv` and `.factor.info.html`: Factor summaries.
    - `.umap.tsv.gz`, `.umap.png`, `.umap.single.prob.png`: Shared UMAP (unless `--skip-umap`).

- Per-sample projection, prefix `samples/<sample>/<sample>.<model_id>`:
    - `.results.tsv.gz`: Factor posteriors of the sample's hexagons.
    - `.umap.tsv.gz`, `.umap.png`, `.umap.single.prob.png`: Per-sample UMAP (unless `--skip-umap`).

- Per-sample decode, prefix `samples/<sample>/<sample>.<model_id>_pixel`:
    - `.bin` and `.index`: Binary pixel-level decoding results and their index.
    - `.pseudobulk.tsv.gz`: Feature-by-factor posterior counts.
    - `.bulk_chisq.tsv`, `.factor.info.tsv`, `.factor.info.html`: Per‑decode DE and factor summaries.
    - `.png`: Image of the pixel-level factors.

- Done flags: `<model_id>.done`, `<model_id>_summary.done`, `<model_id>_each.done`, `<model_id>_pixel_each.done`, `multi_json_each.done`, and per-sample `.done` flags.

- Manifests:
    - `samples/<sample>/ficture.params.json`: Per-sample manifest of the tiled SGE, models, projections, UMAPs and decodes, read by `run_cartload2`. With `--segment-10x`, a `mex` object maps each exported width to its MEX directory.
    - `ficture.multi.params.json`: Shared manifest (see below).

See [`run_ficture2`](./run_ficture2.md#output) for the single‑sample file formats.

---
## Multi-sample manifest

Alongside the per‑sample JSONs, `run_ficture2_multi` writes a top‑level **`ficture.multi.params.json`** that ties the run together. It (a) points to each per‑sample manifest by **relative path**, (b) records the **shared** components — the joint LDA models, shared UMAPs, and multi‑sample hexagon files — and (c) with `--segment-10x`, lists the per‑sample hexagon MEX directories under `mex` (sample → width → directory), so other tools (notably [`run_cartload2_multi`](./run_cartload2_multi.md)) can discover the samples and shared assets without re‑deriving paths. All paths are relative to `--out-dir`, so the directory is self‑contained.

It is written when the Makefile is generated. Entries are merged by `model_id` into an existing manifest, so running once per model (e.g. several projections into the same `--out-dir`) accumulates all models; a `--prepare-only` run overwrites it.

```jsonc
{
  "analysis_type": "multi-sample",
  "n_samples": 2,
  "samples": {
    "s1": "samples/s1/ficture.params.json",
    "s2": "samples/s2/ficture.params.json"
  },
  "shared": {
    "multi_hexagon": { "features": "multi.features.tsv",
                       "hex":  { "12": "multi.hex_12.txt" },
                       "json": { "12": "multi.hex_12.json" } },
    "train_params": [
      { "model_type": "lda", "model_id": "t12_f24", "train_width": 12, "n_factor": 24,
        "cmap": "t12_f24.cmap.tsv", "model_path": "t12_f24.model.tsv",
        "fit_path": "t12_f24.results.tsv.gz",
        "de_path": "t12_f24.bulk_chisq.tsv", "info_path": "t12_f24.factor.info.tsv",
        "umap": { "tsv": "t12_f24.umap.tsv.gz", "png": "t12_f24.umap.png",
                  "ind_png": "t12_f24.umap.single.prob.png" } }
    ]
  },
  "mex": {                                  // only with --segment-10x
    "s1": { "12": "samples/s1/s1.hex_12.mex" },
    "s2": { "12": "samples/s2/s2.hex_12.mex" }
  }
}
```

### Cell analyses

[`run_ficture2_multi_cells`](./run_ficture2_multi_cells.md) (the cell/segmentation‑based analysis) runs on the tiles of a `run_ficture2_multi` output directory and writes the analogous **`ficture.multi.<prefix>.params.json`** for each output prefix (e.g. `ficture.multi.cartloader.params.json`). It points to each `samples/<sample>/ficture.<prefix>.params.json` and records the shared cell components — shared cluster pseudobulk/DE/info, shared heatmap, and shared UMAP/TSNE manifolds — for whichever steps ran.
