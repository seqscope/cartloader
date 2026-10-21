# Multi‑Sample Spatial Asset Packaging

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) calls `run_cartload2_multi` in its `cartload` stage for joint multi-sample runs. This page is for running the packaging step directly, step by step.

## Overview

`run_cartload2_multi` packages a **joint multi‑sample FICTURE run** ([`run_ficture2_multi`](./run_ficture2_multi.md)) into per‑sample PMTiles and catalogs, and writes a **shared multi‑sample catalog**.

It reads the shared manifest [`ficture.multi.params.json`](./run_ficture2_multi.md#multi-sample-manifest) from `--fic-dir` to discover the samples, then fans out one [`run_cartload2`](./run_cartload2.md) per sample from a single Makefile, so samples package **in parallel**. Each sample is written to a self‑contained directory, and a top‑level `multi-catalog.yaml` ties them together.

!!! info "Cell analyses are handled automatically"
    For each sample, every `ficture.<prefix>.params.json` in `samples/<sample_id>/` (e.g. `cartloader`, `xeniumranger`, produced by [`run_ficture2_multi_cells`](./run_ficture2_multi_cells.md)) is auto‑detected and forwarded to that sample's `run_cartload2` as `--in-cell-params`. Shared cell factors come from `ficture.multi.<prefix>.params.json` at the `--fic-dir` root. Background images (e.g. histology) are **not** handled here; add them per sample with [`import_image`](./import_image.md), or use [`run_together`](./run_together.md), which imports images after packaging.

---
## Requirements

- A joint FICTURE output directory (`--fic-dir`) containing **`ficture.multi.params.json`** and per‑sample results under `samples/<sample_id>/`.
- The tools required by [`run_cartload2`](./run_cartload2.md#requirements).

---
## Directory layout

With `--id <multi_id>` (default: the `--out-dir` basename), each sample is packaged into a **self‑contained** directory named `<multi_id>-<sample_id>`:

```
out_dir/
├── multi-catalog.yaml              # links per-sample catalogs + shared factors
├── genes_bin_counts.json           # shared gene->bin assignment (when the manifest lists a joint feature list)
├── t12-f24-model.tsv, t12-f24-rgb.tsv, ...   # shared factor files (see below)
├── t12-f24.shared.done             # one flag per shared factor set
├── <multi_id>-<sample1>/           # self-contained: catalog.yaml + PMTiles
├── <multi_id>-<sample2>/
├── tmp/
└── run_cartload2_multi.mk
```

`multi-catalog.yaml` mirrors the per-sample catalog layout:

- `samples:` — a relative pointer to each `<multi_id>-<sample_id>/catalog.yaml`;
- `assets.sge.counts:` — the shared `genes_bin_counts.json` (present only when the manifest records a joint feature list, `shared.multi_hexagon.features`);
- `assets.factors:` — one entry per shared factor set, keyed by its **hyphenated** id (e.g. `t12-f24`, `cartloader`), model-derived and cell-derived alike, with `post` (the factor×feature matrix: `-model.tsv` for an LDA model, `-pseudobulk.tsv.gz` for cell clusters), `rgb`, `de`, `info`, `umap` (`{pmtiles, png, tsv}`, when available), and `heatmap` (`{pdf, tsv}`, cell factors only). Cell-derived entries include only the files present in the manifest.

The shared files are **materialized into the `out_dir` root** next to `multi-catalog.yaml`, so the catalog references only local basenames.

```yaml
id: <multi_id>
analysis_type: multi-sample
n_samples: 2
samples:
  s1: <multi_id>-s1/catalog.yaml
  s2: <multi_id>-s2/catalog.yaml
assets:
  sge:
    counts: genes_bin_counts.json
  factors:
    t12-f24:                        # model-derived factor set
      post: t12-f24-model.tsv
      rgb:  t12-f24-rgb.tsv
      de:   t12-f24-de.tsv
      info: t12-f24-info.tsv
      umap: { pmtiles: t12-f24-umap.pmtiles, png: t12-f24-umap.png, tsv: t12-f24-umap.tsv.gz }
    cartloader:                     # cell-derived factor set (same keys)
      post: cartloader-pseudobulk.tsv.gz
      rgb:  cartloader-rgb.tsv
      de:   cartloader-de.tsv
      info: cartloader-info.tsv
      umap: { pmtiles: cartloader-umap.pmtiles, png: cartloader-umap.png, tsv: cartloader-umap.tsv.gz }
      heatmap: { pdf: cartloader-heatmap.pdf, tsv: cartloader-heatmap.tsv }
```

Each per-sample `run_cartload2` reuses the shared gene→bin assignment (`--in-bin-json`) and copies the shared UMAP PMTiles from the root (`--reuse-shared-umap-dir`) instead of rebuilding them.

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

Add `--dry-run` to write the Makefile and print the commands without executing them. `multi-catalog.yaml` is written as soon as the command runs, including in dry-run mode.

---
## Parameters

### Input/Output Parameters

- `--fic-dir` (str, **required**): Joint FICTURE directory containing the multi manifest and `samples/<id>/`.
- `--out-dir` (str, **required**): Output root; holds one `<multi_id>-<sample_id>/` per sample plus `multi-catalog.yaml`.
- `--id` (str): Identifier for the run; per‑sample dirs and catalog ids are `<id>-<sample_id>` (default: basename of `--out-dir`).
- `--in-multi-params` (str): Name of the shared manifest under `--fic-dir` (default: `ficture.multi.params.json`).
- `--multi-catalog` (str): Output multi‑catalog filename under `--out-dir` (default: `multi-catalog.yaml`).
- `--out-catalog` (str): Per‑sample catalog filename (default: `catalog.yaml`).

### Run Parameters

- `-j`, `--n-jobs` (int): Parallel jobs (default: `1`). Applied to **both** the outer make (samples in parallel) and each sample's `run_cartload2` (e.g., parallel gene bins), so effective concurrency can reach `n-jobs × n-jobs`; lower it if this oversubscribes the machine.
- `--threads` (int): Threads per job, forwarded to `run_cartload2` (its default: `4`).
- `--dry-run` / `--restart` / `--makefn` (default: `run_cartload2_multi.mk`) / `--log` / `--log-suffix`.

??? note "Auxiliary parameters forwarded to each `run_cartload2`"
    Use defaults unless needed. Forwarded when set: `--in-fic-params`, `--out-fic-assets`, `--rename-x/y`, `--colname-feature/count`, `--out-molecules-id`, `--bin-count`, `--preserve-point-density-thres`, `--sge-scale`, `--use-pmpoint`, `--tile-format-pmpoint`, `--umap-colname-factor/x/y`, `--umap-min/max-zoom`, `--skip-umap`, `--skip-raster`, `--keep-intermediate-files`, and the run options `--restart`, `--threads`, `--log`, `--log-suffix`, `--n-jobs`. (`--max-join-dist-um`, `--join-tile-size`, and `--transparent-below/above` are forwarded but unused by `run_cartload2`.)

    Environment: `--gzip`, `--pmtiles`, `--gdal_translate`, `--gdaladdo`, `--tippecanoe`, `--spatula`, `--pmpoint`, `--ficture2`.

    `--skip-umap` affects the per-sample packaging only; the shared UMAP at the root is still built when the manifest provides one.

    Not forwarded: `--tmp-dir` sets the temporary directory for the root-level shared UMAP only (default: `<out-dir>/tmp`); each sample uses `<sample_dir>/tmp`. Options such as `--max-point-tile-bytes`, `--hex-thres-prob`, or `--use-ficture2-direct-pmtiles` are not exposed here.

    See the [`run_cartload2` reference](./run_cartload2.md) for their meanings and defaults.

---
## Output

Under `--out-dir`:

- `multi-catalog.yaml`: The shared multi‑sample catalog (sample pointers + shared assets).
- `genes_bin_counts.json`: Shared gene→bin assignment and aggregate counts (name follows `--out-molecules-id`; only when the manifest records a joint feature list).
- Shared factor files: `<oid>-model.tsv` or `<oid>-pseudobulk.tsv.gz`, `<oid>-rgb.tsv`, `<oid>-de.tsv`, `<oid>-info.tsv`, `<oid>-umap.{pmtiles,png,tsv.gz}`, and `<oid>-heatmap.{pdf,tsv}` for cell factors; one `<oid>.shared.done` flag per factor set.
- `<multi_id>-<sample_id>/`: One self‑contained directory per sample with its `catalog.yaml`, PMTiles, and asset JSONs (see [`run_cartload2` → Output](./run_cartload2.md#output)).
- `tmp/` and `run_cartload2_multi.mk`.

For the upstream pipeline, see [`run_ficture2_multi`](./run_ficture2_multi.md); for a fully orchestrated end‑to‑end run, see [`run_together`](./run_together.md).
