# Spatial Asset Packaging

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) calls `run_cartload2` (single sample) or [`run_cartload2_multi`](./run_cartload2_multi.md) (joint runs) in its `cartload` stage. This page is for running the packaging step directly, step by step.

## Overview

`run_cartload2` packages SGE data and spatial factor output from [FICTURE analysis](./run_ficture2.md) into standardized, spatially indexed, and storage‑efficient [PMTiles](https://github.com/protomaps/PMTiles), a web‑native tiling format, and writes a `catalog.yaml` describing every layer. These outputs are optimized for interactive web visualization (e.g., in CartoScope) and data sharing.

---
## Requirements

- A completed FICTURE run (`--fic-dir`) from [`run_ficture2`](./run_ficture2.md) (or one sample directory of [`run_ficture2_multi`](./run_ficture2_multi.md)) with its `ficture.params.json`. Its `in_sge` block supplies the tiled transcript (`<in_tiled>.tsv` + `<in_tiled>.index`), feature, and minmax files.
- (Optional) Cell-analysis parameter JSONs from [`run_ficture2_multi_cells`](./run_ficture2_multi_cells.md) (`--in-cell-params`).
- (Optional) Asset JSON/YAML files from [`import_image`](./import_image.md), [`import_*_cell`](./import_cell.md), or [`import_visiumhd_square`](./import_square.md).
- Tools (defaults point to the bundled submodules):
    - always: `spatula`, `tippecanoe`, and punkst (`--ficture2`; `<ficture2>/bin/punkst` must exist even without FICTURE results);
    - unless `--skip-raster`: `pmtiles`;
    - with `--use-image-png2pmtiles`: `gdal_translate`, `gdaladdo`;
    - with `--use-pmpoint` or `--use-ficture2-direct-pmtiles`: `pmpoint`.

!!! warning "Packaging without FICTURE (`--sge-dir`)"
    With `--sge-dir`, the SGE assets JSON must provide `transcript`, `feature`, and `minmax`, and `<transcript>.tsv` and `<transcript>.index` (a punkst tiled transcript) must exist. The `sge_assets.json` written by [`sge_convert`](./sge_convert.md) points to `transcripts.unsorted.tsv.gz`, which does not meet this requirement, so `--sge-dir` alone currently fails on standard `sge_convert` output. If both `--sge-dir` and `--fic-dir` are given, `--sge-dir` is used when its JSON has all three keys. Use `--fic-dir` (a `run_ficture2` output) as the SGE source.

---
## Actions

!!! warning "Action Specifications"
    All applicable steps run by default; which ones apply depends on the [inputs](#inputoutput-parameters) you provide.

### SGE Packaging
Rasterizes the tiled transcripts into light and dark mono PMTiles (`sge-mono-light.pmtiles`, `sge-mono-dark.pmtiles`), unless `--skip-raster`.

### FICTURE Integration (`--fic-dir`)
For each trained model in `ficture.params.json`:

- Converts the hexagon-level factor results (`<model_id>.results.tsv.gz`) into a vector PMTiles layer (`<model-id>.pmtiles`).
- Copies the model matrix, color table, DE, and factor info under renamed files (see [Output](#output)).
- Converts UMAP embeddings into PMTiles (unless `--skip-umap`).
- Converts each decoded pixel image into a raster PMTiles layer (`<decode-id>-pixel-raster.pmtiles`), unless `--skip-raster`.
- Annotates every molecule with the decoded factors of its pixel (`punkst tile-op`), then splits the joined molecules into gene-bin point PMTiles.

Without decoded factors (e.g., a `run_ficture2_multi --prepare-only` manifest), the tiled transcripts are converted to point PMTiles directly.

### Cell Analysis Integration (`--in-cell-params`)
For each cell analysis: builds cell-centroid and boundary PMTiles, copies its DE, info, pseudobulk, and heatmap files, converts its pixel image into a raster layer, and joins its decoded pixels into the molecule layers.

### Additional Assets (`--cell-assets`, `--background-assets`, `--square-assets`)
Adds the provided assets to `catalog.yaml`; referenced files are copied into `--out-dir` when they live elsewhere.

### Catalog and Metadata
Writes `ficture_assets.json` (with `--fic-dir`) and `catalog.yaml`. An existing `catalog.yaml` in `--out-dir` is updated in place (assets added or replaced) rather than overwritten.

---
## Example Usage

```bash
DATA_ID="dataset-id"               ## replace dataset-id with the id for your dataset

cartloader run_cartload2 \
    --fic-dir /path/to/run_ficture2/results \
    --out-dir /path/to/output/directory \
    --id ${DATA_ID} \
    --colname-count count \
    --n-jobs 20  \
    --threads 20
```

Tool paths default to the bundled submodules; override them (e.g., `--spatula`, `--tippecanoe`, `--pmtiles`, `--ficture2`) only if your binaries live elsewhere.

---
## Parameters

Below are the core parameters. See more details in the collapsible "Auxiliary Parameters" section.

### Input/Output Parameters
Provide the SGE through `--fic-dir` (recommended) or `--sge-dir` (see the warning above).

* `--out-dir` (str, required): Output directory for PMTiles, assets JSON, and catalog YAML.
* `--fic-dir` (str): `run_ficture2` output directory; must include the `--in-fic-params` JSON/YAML.
* `--sge-dir` (str): `sge_convert` output directory; must include the `--in-sge-assets` JSON/YAML.
* `--in-cell-params` (list of str): Cell-analysis parameter JSON/YAML files from `run_ficture2_multi_cells` (each must contain a `cell_params` block).
* `--cell-assets` (list of str): Cell asset JSON/YAML files from `import_xenium_cell` / `import_visiumhd_cell`.
* `--background-assets` (list of str): Background assets from `import_image`, as JSON/YAML files or inline `id:path` / `id1:id2:path` specs.
* `--square-assets` (list of str): Square asset JSON/YAML files from `import_visiumhd_square`.

### Dataset ID and Descriptions

* `--id` (str, required): Dataset identifier (no whitespace; use `-`).
* `--title` (str): Title for the catalog (default: `--id`).
* `--desc` (str): Short description for the catalog.

??? note "Auxiliary Parameters"
    Recommend to use the default values; override only if needed.

    **Input/Output Names**:

    * `--in-sge-assets` (str): SGE assets JSON/YAML under `--sge-dir` (default: `sge_assets.json`).
    * `--in-fic-params` (str): FICTURE parameter JSON/YAML under `--fic-dir` (default: `ficture.params.json`).
    * `--out-fic-assets` (str): Output FICTURE assets JSON/YAML under `--out-dir` (default: `ficture_assets.json`).
    * `--out-catalog` (str): Output catalog YAML under `--out-dir` (default: `catalog.yaml`).
    * `--out-molecules-id` (str): Base name of the molecule PMTiles files (default: `genes`).
    * `--tmp-dir` (str): Temporary directory (default: `<out-dir>/tmp`).
    * `--keep-intermediate-files` (flag): Keep intermediate files.

    **Molecule (Point) PMTiles**:

    * `--colname-feature` (str): Feature column name (default: `gene`).
    * `--colname-count` (str): Count column name (default: `count`).
    * `--rename-x` / `--rename-y` (str): Column renaming rules for tippecanoe (default: `x:lon` / `y:lat`).
    * `--bin-count` (int): Number of gene bins for the point PMTiles (default: `50`).
    * `--max-point-tile-bytes` (int): Maximum tile size in bytes for point layers (default: `10000000`).
    * `--max-point-feature-counts` (int): Maximum features per tile for point layers (default: `1000000`).
    * `--preserve-point-density-thres` (int): Tippecanoe point-density preservation threshold (default: `1024`).
    * `--replace-features` (str): Shared feature list (e.g., `multi.features.tsv`); a header-normalized copy (`multi.features.rehdr.tsv`) is used for gene-to-bin assignment so batches sharing the list bin genes identically.
    * `--in-bin-json` (str): Precomputed gene→bin assignment JSON (a `_bin_counts.json`); reused instead of deriving bins, with counts refreshed from this sample. Takes precedence over `--replace-features`. Supplied automatically by `run_cartload2_multi`.
    * `--use-pmpoint` (flag): Build point PMTiles with `pmpoint` instead of tippecanoe.
    * `--tile-format-pmpoint` (str): `MVT` or `MLT` with `--use-pmpoint` (default: `MVT`).
    * `--pmpoint-compression-scale` (float): Extra pmpoint compression scale (default: `10.0`).
    * `--use-ficture2-direct-pmtiles` (flag): Package joined molecules directly with `punkst tile-op` + `pmpoint` instead of `run_tsv2pmtiles`; applies only when decoded factors are joined, and requires a positive `--bin-count`.
    * `--point-min-zoom` / `--point-max-zoom` (int): Zoom range for `--use-ficture2-direct-pmtiles` only (default: `10` / `18`).

    **Polygon PMTiles (hexagons, cells, boundaries)**:

    * `--max-polygon-tile-bytes` (int): Maximum tile size in bytes (default: `50000000`).
    * `--max-polygon-feature-counts` (int): Maximum features per tile (default: `5000000`).
    * `--hex-thres-prob` (float): Minimum per-factor probability stored in hexagon PMTiles (default: `0.0001`).

    **Raster PMTiles**:

    * `--skip-raster` (flag): Skip SGE and factor raster layers (no `pmtiles` binary required).
    * `--sge-scale` (int): Units per pixel for rasterized images (default: `1`).
    * `--raster-min-zoom` / `--raster-max-zoom` (int): Zoom range for decoded-factor raster layers (default: `6` / `18`).
    * `--use-image-png2pmtiles` (flag): Use the GDAL-based `image_png2pmtiles` path instead of `punkst image2pmtiles`.

    **UMAP**:

    * `--skip-umap` (flag): Skip UMAP PMTiles and copied UMAP files even if UMAP inputs exist.
    * `--umap-colname-factor` (str): Dominant-factor column in the UMAP TSV (default: `topK`).
    * `--umap-colname-x` / `--umap-colname-y` (str): UMAP coordinate columns (default: `UMAP1` / `UMAP2`).
    * `--umap-min-zoom` / `--umap-max-zoom` (int): UMAP PMTiles zoom range (default: `0` / `18`).
    * `--reuse-shared-umap-dir` (str): Directory with prebuilt shared `<factor>-umap.pmtiles` to copy instead of re-running tippecanoe (set by `run_cartload2_multi`).

    **Accepted but unused**: `--max-join-dist-um`, `--join-tile-size`, `--transparent-below`, `--transparent-above` are parsed but have no effect.

    **Environment Parameters** (defaults are submodule paths unless noted):

    * `--spatula` (str): `spatula` binary (default: `<cartloader_dir>/submodules/spatula/bin/spatula`).
    * `--tippecanoe` (str): `tippecanoe` binary (default: `<cartloader_dir>/submodules/tippecanoe/tippecanoe`).
    * `--pmtiles` (str): go-pmtiles binary (default: `<cartloader_dir>/submodules/pmtiles/pmtiles`).
    * `--ficture2` (str): punkst repository (default: `<cartloader_dir>/submodules/punkst`).
    * `--pmpoint` (str): `pmpoint` binary (default: `<cartloader_dir>/submodules/pmpoint/bin/pmpoint`).
    * `--gdal_translate` / `--gdaladdo` (str): GDAL binaries on `PATH` (default: `gdal_translate` / `gdaladdo`).
    * `--gzip` (str): gzip-compatible binary on `PATH`, e.g., `pigz` (default: `gzip`).

    **Run Parameters**:

    * `--dry-run` (flag): Write the Makefile and print the commands without executing them.
    * `--restart` (flag): Ignore existing outputs and rerun all steps.
    * `--makefn` (str): Makefile name under `--out-dir` (default: `run_cartload2.mk`).
    * `--n-jobs` (int): Number of parallel jobs (default: `1`).
    * `--threads` (int): Max threads per job (default: `4`).
    * `--log` (flag): Write a log to `<out-dir>_cartload<log-suffix>` (next to `--out-dir`).
    * `--log-suffix` (str): Log filename suffix (default: `.log`).

---
## Output

All files are written under `--out-dir`. FICTURE ids are hyphenated in file names (e.g., model `t12_f24` → `t12-f24`, decode `t12_f24_pixel` → `t12-f24-pixel`).

### Rasterized Transcript-level SGE
* `sge-mono-dark.pmtiles`, `sge-mono-light.pmtiles`: Rasterized SGE for dark and light backgrounds.

### Spatial Factor Layers (per model)
* `<model-id>.pmtiles`: Hexagon-level factor probabilities (vector tiles).
* `<model-id>-model.tsv`, `<model-id>-rgb.tsv`, `<model-id>-bulk-de.tsv`, `<model-id>-info.tsv`: Copied model matrix, color table, DE, and factor info (`<model-id>-factor-map.tsv` when present).
* `<model-id>-umap.pmtiles`, `<model-id>-umap.tsv.gz`, `<model-id>.umap.png`, `<model-id>.umap.single.prob.png`: UMAP layer and copies (multi-sample models also get `<model-id>-shared-umap.*` and `<model-id>-shared.umap*.png`).

### Decoded Factor Layers (per decode)
* `<decode-id>-pixel-raster.pmtiles`: Rasterized pixel-level factor map.
* `<decode-id>-bulk-de.tsv`, `<decode-id>-pseudobulk.tsv.gz`, `<decode-id>-info.tsv`: Copied decode DE, pseudobulk, and factor info.

### Cell Analysis Layers (with `--in-cell-params`)
* `<cell-id>-cells.pmtiles`, `<cell-id>-boundaries.pmtiles`, `<cell-id>-pixel-raster.pmtiles`: Cell centroids, boundaries, and pixel raster.
* `<cell-id>-rgb.tsv`, `<cell-id>-bulk-de.tsv`, `<cell-id>-info.tsv`, `<cell-id>-pseudobulk.tsv.gz`, `<cell-id>-heatmap.pdf`, `<cell-id>-heatmap.tsv` (plus `-shared-*` copies for multi-sample analyses) and UMAP files as above.

### Molecule PMTiles
* `transcripts_pixel_joined.tsv`: Molecules annotated with decoded factors (uncompressed; kept after the run).
* `genes_bin*.pmtiles`, `genes_all.pmtiles`: Gene-bin point layers and the all-gene layer (`genes_all.pmtiles` only when factors are joined).
* `genes_index.tsv`, `genes_pmtiles_index.tsv`, `genes_bin_counts.json`: Bin index, PMTiles index, and per-bin gene counts (prefix set by `--out-molecules-id`).

### Summary Files
* `ficture_assets.json`: FICTURE assets per trained model.
* `catalog.yaml`: Catalog of all layers for CartoScope.
* `run_cartload2.mk` and `*.done` flag files (e.g., `<model-id>.done`, `<decode-id>.done`, `ficture.done`, `sge-mono.done`).
