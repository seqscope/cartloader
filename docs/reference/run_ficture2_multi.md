# Multi-Sample Spatial Factor Inference Analysis with FICTURE

## Overview

`run_ficture2_multi` orchestrates joint FICTURE analysis on >=2 samples. It tiles each sample, builds joint hex grids (one or more widths), trains LDA models cross samples per width/factor count, decodes per sample, and writes per‑sample JSON manifests summarizing results.

Use this when processing multiple samples together to learn shared spatial factors and generate per‑sample outputs in a single, parallelizable pipeline.

---
## Requirements

- TSV file where the first two columns per line are the sample ID and the path to the transcript‑indexed SGE file  
- FICTURE2 repository with `bin/punkst` and Python utilities (e.g., `ext/py/factor_report.py`)
- Pre-installed tools: `gzip`, `python`, `punkst`, `spatula`
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
## Actions

All actions run by default (UMAPs can be skipped with `--skip-umap`):

- Multisample prepare: tiles inputs, builds joint hex grids at requested `--width` values
- LDA training: trains LDA models for each `(width, n-factor)` pair
- Decode: applies trained models per sample; produces pixel‑level factors and summaries
- Write per‑sample JSON: consolidates paths and metadata for downstream consumption

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

### Training / Decoding

- `--n-factor` (str): Comma‑separated factor counts for training (e.g., `12,24`).
- `--anchor-res` (int): Anchor resolution used in decode IDs (see outputs).
- `--cmap-file` (str, defaults to [fixed_color_map_256.tsv](../../assets/fixed_color_map_256.tsv)): Colormap TSV used to colorize factors.
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

- Per‑sample JSONs
    - `samples/<sample>/ficture.params.json`: Consolidates sample feature paths, LDA and decode outputs for downstream steps.

---
## See Also

- Reference: `run_ficture2.md` — Single‑sample FICTURE2 runner and file formats
