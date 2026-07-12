# Multi‑Sample Spatial Asset Packaging

## Overview

`run_cartload2_multi` packages a **joint multi‑sample FICTURE run** ([`run_ficture2_multi`](./run_ficture2_multi.md)) into per‑sample PMTiles and catalogs, and writes a **shared multi‑sample catalog**.

It reads the shared manifest [`ficture.multi.params.json`](./run_ficture2_multi.md#multi-sample-manifest) from `--fic-dir` to discover the samples, then fans out one [`run_cartload2`](./run_cartload2.md) per sample as a single Makefile so they package **in parallel** and fault‑tolerantly. Each sample is written to a self‑contained directory, and a top‑level `multi-catalog.yaml` ties them together.

!!! info "Cell analyses are handled automatically"
    For each sample, any `ficture.<prefix>.params.json` present (e.g. `cartloader`, `xeniumranger`, produced by [`run_ficture2_multi_cells`](./run_ficture2_multi.md)) is auto‑detected and forwarded to that sample's `run_cartload2` as `--in-cell-params`. Background images (e.g. histology) are **not** handled here — add them per sample with `import_image`, or use [`run_together`](./run_together.md), which runs image import after packaging.

---
## Requirements

- A joint FICTURE output directory (`--fic-dir`) containing **`ficture.multi.params.json`** and per‑sample results under `samples/<sample_id>/`.
- CLI tools used by CartLoad2: `tippecanoe`, `gdal`, `pmtiles`, `gzip`/`pigz`, `spatula` (and `pmpoint` with `--use-pmpoint`).

---
## Directory layout

With `--id <multi_id>` (default: the `--out-dir` basename), each sample is packaged into a **self‑contained** directory named `<multi_id>-<sample_id>`, so the whole `--out-dir` uploads to S3 as one unit:

```
out_dir/
├── multi-catalog.yaml              # links per-sample catalogs + shared factors (copied to root)
├── <multi_id>-<sample1>/           # self-contained: catalog.yaml + PMTiles
├── <multi_id>-<sample2>/
└── run_cartload2_multi.mk
```

The `multi-catalog.yaml` records:

- `samples:` — a relative pointer to each `<multi_id>-<sample_id>/catalog.yaml` (the **only** pointers into sample subdirectories);
- `factors:` — one entry per shared factor set, **model-derived and cell-derived alike** (e.g. `t12_f24`, `cartloader`), all with the same keys: `post` (factor×feature matrix — the LDA `model.tsv` or the cell pseudobulk), `rgb`, `de`, `info`, `umap` (`{pmtiles, png, tsv}`), and `heatmap` (`{pdf, tsv}`, cell factors only).

Those shared files are **copied into the `out_dir` root** next to `multi-catalog.yaml` (materialized once from sample‑1's already-processed outputs), so the catalog references only local basenames and `out_dir` is fully self‑contained.

```yaml
id: <multi_id>
analysis_type: multi-sample
n_samples: 2
samples:
  s1: <multi_id>-s1/catalog.yaml
  s2: <multi_id>-s2/catalog.yaml
factors:
  t12_f24:                        # model-derived factor set
    post: t12-f24-post.tsv
    rgb:  t12-f24-rgb.tsv
    de:   t12-f24-de.tsv
    info: t12-f24-info.tsv
    umap: { pmtiles: t12-f24-umap.pmtiles, png: t12-f24-umap.png, tsv: t12-f24-umap.tsv.gz }
  cartloader:                     # cell-derived factor set (same keys)
    post: cartloader-post.tsv.gz
    rgb:  cartloader-rgb.tsv
    de:   cartloader-de.tsv
    info: cartloader-info.tsv
    umap: { pmtiles: cartloader-umap.pmtiles, png: cartloader-umap.png, tsv: cartloader-umap.tsv.gz }
    heatmap: { pdf: cartloader-heatmap.pdf, tsv: cartloader-heatmap.tsv }
```

---
## Example Usage

```bash
cartloader run_cartload2_multi \
  --fic-dir /path/to/run_ficture2_multi/out \
  --out-dir /path/to/out/mydataset \
  --id      mydataset \
  --use-pmpoint --bin-count 500 \
  -j 10 --threads 24
```

Add `--dry-run` to write the Makefile and print the per‑sample commands without executing.

---
## Parameters

### Input/Output Parameters

- `--fic-dir` (str, **required**): Joint FICTURE directory containing the multi manifest and `samples/<id>/`.
- `--out-dir` (str, **required**): Output root; holds one `<multi_id>-<sample_id>/` per sample plus `multi-catalog.yaml`.
- `--id` (str): Identifier for the run; per‑sample dirs are `<id>-<sample_id>` (default: basename of `--out-dir`).
- `--in-multi-params` (str): Name of the shared manifest under `--fic-dir` (default: `ficture.multi.params.json`).
- `--multi-catalog` (str): Output multi‑catalog filename under `--out-dir` (default: `multi-catalog.yaml`).
- `--out-catalog` (str): Per‑sample catalog filename (default: `catalog.yaml`).

### Run Parameters

- `-j`, `--n-jobs` (int): Number of samples to package **in parallel** (default: 1).
- `--threads` (int): Threads per job, forwarded to `run_cartload2`.
- `--dry-run` / `--restart` / `--makefn` / `--log` / `--log-suffix`.

??? note "Auxiliary parameters forwarded to each `run_cartload2`"
    Use defaults unless needed. Packaging: `--in-fic-params`, `--out-fic-assets`, `--rename-x/y`, `--colname-feature/count`, `--out-molecules-id`, `--max-join-dist-um`, `--join-tile-size`, `--bin-count`, `--preserve-point-density-thres`, `--sge-scale`, `--use-pmpoint`, `--tile-format-pmpoint`, `--umap-colname-*`, `--umap-min/max-zoom`, `--skip-umap`, `--skip-raster`, `--tmp-dir`, `--keep-intermediate-files`, `--transparent-below/above`.

    Environment: `--gzip`, `--pmtiles`, `--gdal_translate`, `--gdaladdo`, `--tippecanoe`, `--spatula`, `--pmpoint`, `--ficture2`.

    See the [`run_cartload2` reference](./run_cartload2.md) for their meanings and defaults.

---
## Output

Under `--out-dir`:

- `multi-catalog.yaml` — the shared multi‑sample catalog (pointers + shared assets).
- `<multi_id>-<sample_id>/catalog.yaml` + PMTiles and asset JSONs, one self‑contained directory per sample (see [`run_cartload2` → Output](./run_cartload2.md#output)).
- `run_cartload2_multi.mk` — the Makefile capturing all per‑sample targets.

For the upstream pipeline, see [`run_ficture2_multi`](./run_ficture2_multi.md); for a fully orchestrated end‑to‑end run, see [`run_together`](./run_together.md).
