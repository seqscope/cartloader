# Multi-Sample Spatial Factor Inference Analysis with FICTURE

## Overview

`run_ficture2_multi` orchestrates joint FICTURE analysis on >=2 samples. It tiles each sample, builds joint hex grids (one or more widths), trains LDA models cross samples per width/factor count, decodes per sample, and writes per‑sample JSON manifests summarizing results.

Use this when processing multiple samples together to learn shared spatial factors and generate per‑sample outputs in a single, parallelizable pipeline.

!!! question
    - [What is `FICTURE` and `punkst`?](../faq/ficture.md)
    - [Why use Multi‑Sample FICTURE Analysis?](../faq/why_multi-sample.md)

---
## Requirements

- TSV file where the first two columns per line are the sample ID and the path to the transcript‑indexed SGE file  
- FICTURE2 repository with `bin/punkst` and Python utilities (e.g., `ext/py/factor_report.py`)
- Pre-installed tools: `gzip`, `python`, `punkst`, `spatula`
---
## Actions

All actions run by default (UMAPs can be skipped with `--skip-umap`):

- Multisample prepare: tiles inputs, builds joint hex grids at requested `--width` values
- Hexagon MEX export (optional, `--segment-10x`): converts each per-sample hexagon file into a 10x MEX directory
- LDA training: trains LDA models for each `(width, n-factor)` pair
- Decode: applies trained models per sample; produces pixel‑level factors and summaries
- Write per‑sample JSON: consolidates paths and metadata for downstream consumption

`--prepare-only` stops after the first action: it tiles the inputs and builds the hex grids, then writes each sample's `ficture.params.json` with only `in_sge` (tiled prefix, feature list, coordinate range) and an empty `train_params`. No model is trained, projected, or decoded, and no UMAP is built. This is what lets a dataset be **packaged without any factor analysis** — `run_cartload2` reads the tiled TSV as its molecule source, so it produces transcripts + raster + images and a catalog with no factor layers. Re-running without `--prepare-only` in the same `--out-dir` adds models later and re-uses the existing tiles/hexagons. See [`run_together --no-ficture`](./run_together_inputs.md#ficture-mode-de-novo-vs-projection).

`--segment-10x` additionally exports every per-sample hexagon file as a **10x MEX** directory, `samples/<sample>/<sample>.hex_<width>.mex/` (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`), via `spatula sptsv2mex`. The MEX is built from the very same `samples/<sample>/<sample>.hex_<width>.txt`/`.json` that the factor analysis uses — the hexagons are neither re-generated nor re-filtered, so the MEX contains exactly the hexagons that survive `--min-ct-per-unit-hexagon`. Barcodes are the hexagon centers as `x:y` (µm); features carry the gene name as both id and name, typed `Gene Expression`. By default every `--width` is exported; `--segment-width-10x` picks the widths instead (a width not in `--width` is simply added to the hexagon step, so its hexagon files get built too). The export depends only on the prepare step, so it also works with `--prepare-only` — the way to get hexagon MEX files without any factor analysis. Each exported directory is recorded under `mex` (keyed by width) in the per-sample `ficture.params.json` and in `ficture.multi.params.json` (keyed by sample, then width).

---
## Example Usage

```bash
cartloader run_ficture2_multi \
  --in-list /path/to/samples.tsv \
  --out-dir /path/to/out/ficture_multi \
  --width 12 \
  --n-factor 12,24 \
  --exclude-feature-regex "^(Blank-.*$)" \
  --redo-merge-units \
  --spatula /path/to/spatula/binary \
  --ficture2 /path/to/FICTURE2 \
  --cmap-file /path/to/colormap.tsv \
  --threads 8 \
  --n-jobs 4
```

---
## Parameters

Below are the core arguments you’ll typically set. Flag names and behavior follow `run_ficture2_multi.py`.

### Input/Output

- `--in-list` (str, required): TSV with at least two columns per line: `<sample_id>` and `<path_to_transcript_TSV>` .
- `--out-dir` (str, required): Output directory.
- `--out-json` (str): Top‑level JSON (defaults to `<out-dir>/ficture.params.json`).

### Hex Tiling / Preprocessing

- `--colidx-x`, `--colidx-y` (int, 1‑based): X/Y column indices in the transcript TSV.
- `--colidx-feature`, `--colidx-count` (int, 1‑based): Feature and count indices.
- `--tile-size` (int): Tile size for preprocessing.
- `--tile-buffer` (int): Buffer size for tiling.
- `--width` (str): Comma‑separated hex widths in µm (e.g., `8,16`).
- `--min-count` (int): Minimum count per unit hexagon.
- `--min-total-count-per-sample` (int): Minimum per‑sample transcript count to retain in the joint set.
- `--include-feature-regex` / `--exclude-feature-regex` (str): Feature filters.
- `--redo-merge-units` (flag): Rebuild merged units per width (temporary bug workaround).
- `--prepare-only` (flag): Tile and build hexagons only; write manifests with no models (`--n-factor` not required, UMAP forced off).
- `--segment-10x` (flag): Also export each per-sample hexagon file as a 10x MEX directory (`samples/<sample>/<sample>.hex_<width>.mex/`) with `spatula sptsv2mex`; same hexagons and `--min-ct-per-unit-hexagon` filter as the analysis. Works with `--prepare-only`.
- `--segment-width-10x` (str): Comma‑separated widths (µm) to export with `--segment-10x` (default: the `--width` list). Extra widths are added to the hexagon step.

### Training / Decoding

- `--n-factor` (str): Comma‑separated factor counts for training (e.g., `12,24`).
- `--anchor-res` (int): Anchor resolution used in decode IDs (see outputs).
- `--cmap-file` (str, defaults to [fixed_color_map_256.tsv](https://github.com/seqscope/cartloader/blob/main/assets/fixed_color_map_256.tsv): Colormap TSV used to colorize factors.
- `--umap` (flag): Generate UMAP embeddings/plots for each LDA model (on by default).
- `--skip-umap` (flag): Skip UMAP generation (overrides `--umap`).

### Run Options

- `--dry-run` (flag): Generate Makefile and print commands only.
- `--restart` (flag): Ignore existing outputs and rerun steps.
- `--threads` (int): Max threads per job (default: 8).
- `--n-jobs` (int): Parallel jobs for the Makefile.

### Environment / Tools

- `--ficture2` (str, required): Path to FICTURE2 repo containing `bin/punkst`.
- `--python` (str): Python executable (used for reporting utilities).
- `--gzip` (str): Path to `gzip` binary.

---
## Output

Outputs are written under `--out-dir`.

- Multisample prepare (joint):
    - `multi.features.tsv`: Joint (tiled) features table.
    - `multi.hex_<width>.txt` and `multi.hex_<width>.json`: Hex grid and metadata per width.

- LDA training (per width × n‑factor): prefix `t{width}_f{n_factor}`
    - `.model.tsv`: Topic–feature weights (factors).
    - `.results.tsv.gz`: Posterior per unit (hex/pixel) with top factor columns.
    - `.bulk_chisq.tsv`: Per‑factor differential feature table.
    - `.factor.info.tsv` and optional `.factor.info.html`: Factor summaries and colors.
    - `.umap.tsv.gz`, `.umap.png`, `.umap.single.prob.png`: UMAP coordinates and plots for factors (written unless `--skip-umap`).

- Decode (per sample × model/width): prefix `<sample>.<decode_id>`
    - `.tsv.gz`: Pixel‑level decode with posterior/assignments.
    - `.png`: Quick‑look image per decode.
    - `.pseudobulk.tsv.gz`: Aggregated counts by factor.
    - `.bulk_chisq.tsv`, `.factor.info.tsv`: Per‑decode summaries.

- Hexagon MEX export (`--segment-10x`; per sample × width)
    - `samples/<sample>/<sample>.hex_<width>.mex/`: `barcodes.tsv.gz` (hexagon centers as `x:y`), `features.tsv.gz`, `matrix.mtx.gz` — the 10x MEX form of `samples/<sample>/<sample>.hex_<width>.txt`.

- Per‑sample JSONs
    - `samples/<sample>/ficture.params.json`: Consolidates sample feature paths, LDA and decode outputs for downstream steps. With `--segment-10x`, a `mex` object maps each exported width to its MEX directory.

- Shared multi‑sample manifest
    - `ficture.multi.params.json`: see below.

See `run_ficture2.md` for single‑sample formats and file details.

---
## Multi-sample manifest

Alongside the per‑sample JSONs, `run_ficture2_multi` writes a top‑level **`ficture.multi.params.json`** that ties the run together. It (a) points to each per‑sample manifest by **relative path**, (b) records the **shared** components — the joint LDA models, shared UMAPs, and multi‑sample hexagon files — and (c) with `--segment-10x`, lists the per‑sample hexagon MEX directories under `mex` (sample → width → directory), so other tools (notably [`run_cartload2_multi`](./run_cartload2_multi.md)) can discover the samples and shared assets without re‑deriving paths. All paths are relative to `--out-dir`, so the directory is self‑contained.

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

`run_ficture2_multi_cells` (the cell/segmentation‑based decode) writes the analogous **`ficture.multi.<prefix>.params.json`** for each output prefix (e.g. `ficture.multi.cartloader.params.json`). It points to each `samples/<sample>/ficture.<prefix>.params.json` and records the shared cell components — shared cluster pseudobulk/DE/info, shared heatmap, and shared UMAP/TSNE manifolds — for whichever steps ran.
