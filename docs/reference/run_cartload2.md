# Spatial Asset Packaging

## Overview

Following spatial factor inference via FICTURE analysis, `CartLoader` provides the `run_cartload2` module to package SGE data and spatial factor output from [FICTURE analysis](./run_ficture2.md) into standardized, spatially indexed, and storage‑efficient [PMTiles](https://github.com/protomaps/PMTiles), a web‑native tiling format. These PMTiles outputs are optimized for downstream analysis, interactive web visualization (e.g., in CartoScope), and data sharing across platforms.

---
## Requirements
- An SGE in the unified format (from [SGE format conversion](./sge_convert.md)) with a metadata file describing paths to SGE assets (`sge_assets.json`)
- (Optional) A completed FICTURE run from [`run_ficture2`](./run_ficture2.md), including pixel-level decoding outputs and a metadata file describing the input-output structure (`ficture.params.json`)
- (Optional) One or more deployed images with corresponding metadata file(s)
- (Optional) A set of deployed cell-based analysis results and its metadata file
- Pre-installed CLI tools: `tippecanoe`, `gdal_translate`, `gdaladdo`, `pmtiles`, `spatula`, `gzip`.

---
## Example Usage

```bash
DATA_ID="dataset_id"               ## replace dataset_id with the id for your dataset

cartloader run_cartload2 \
    --makefn run_cartload2.mk \
    --fic-dir /path/to/run_ficture2/results \
    --out-dir /path/to/output/directory \
    --id ${DATA_ID} \
    --colname-count count \
    --n-jobs 20  \
    --threads 20 \
    --spatula /path/to/spatula/binary \
    --pmtiles /path/to/pmtiles/binary \
    --tippecanoe /path/to/tippecanoe/binary
```

---
## Actions

!!! warning "Action Specifications"
    SGE packaging and catalog preparation run by default. Enable optional integration actions using [input/output parameters](#inputoutput-parameters).

### SGE Packaging
Converts transcript‑level SGE to raster PMTiles, including light and dark modes.

### FICTURE Integration (`--fic-dir`)
If FICTURE results is provided via `--fic-dir`:
  - Converts topic proportions (`*.results.tsv.gz`) into vector PMTiles for spatial factors.
  - Generates raster overlays from decoded pixel‑level outputs.
  - Builds a joined molecule–factor matrix by associating decoded pixels with molecules based on spatial proximity, then converts it into multi‑feature PMTiles.

### Additional Assets Integration (`--cell-assets` or `--background-assets`)
Includes provided cell assets and background/basemap PMTiles in `catalog.yaml`; copies asset files into the output directory when they reside elsewhere.

### Catalog and Metadata Preparation
Writes a FICTURE assets JSON (when FICTURE integration is set) and a final `catalog.yaml` listing all layers.

---
## Parameters

Below are the core parameters. See more details in the collapsible "Auxiliary Paramaters" section.

### Input/Output Parameters
Must use `--sge-dir` or `--fic-dir` to provide SGE data.

* `--sge-dir` (str): Path to the input Directory from `sge_convert`; must include an SGE assets JSON (`--in-sge-assets`). 
* `--fic-dir` (str): Path to the input Directory from `run_ficture2`; must include `--in-fic-params` JSON/YAML.
* `--cell-assets` (list of str): Optional cell asset JSON/YAML files (e.g., from `import_*_cell`) to include.
* `--background-assets` (list of str): Optional basemap assets (e.g., from cell `import_image`); JSON/YAML or inline specs `id:path` or `id1:id2:path`.
* `--square-assets` (list of str): Optional square asset JSON/YAML files (e.g., from `import_visiumhd_square`) to include.
* `--out-dir` (str): Path to the output directory for PMTiles, assets JSON, and catalog YAML.

### Dataset ID and Descriptions

* `--id` (str): Unique dataset identifier (avoid whitespace; use `-`).
* `--title` (str): Optional title for the output assets.
* `--desc` (str): Optional short description for the output assets.

??? note "Auxiliary Parameters"
    Recommend to use the default values; override only if needed.

    **Auxiliary Conversion Parameters**:
    
    * `--in-sge-assets` (str): File name of SGE assets JSON/YAML under `--sge-dir` specifying paths to transcript, feature, and minmax files (default: `sge_assets.json`).
    * `--in-fic-params` (str): File name of input JSON/YAML with SGE paths and FICTURE parameters under `--fic-dir`(default: `ficture.params.json`).
    * `--out-fic-assets` (str): File name of output JSON/YAML file to write FICTURE assets (default: `ficture_assets.json`).
    * `--out-catalog` (str): File name of output YAML file for assets (default: `catalog.yaml`).
    * `--rename-x` (str): Column renaming rule for X axis in `tippecanoe` (Default: x:lon).
    * `--rename-y` (str): Column renaming rule for Y axis in `tippecanoe` (Default: y:lat).
    * `--colname-feature` (str): Column name for gene/feature name (default: `gene`).
    * `--colname-count` (str): Column name for feature count (default: `count`).
    * `--out-molecules-id` (str): Base name for output molecule PMTiles files (default: `genes`).
    * `--max-join-dist-um` (float): Maximum join distance (µm) between molecules and pixels (default: `0.1`).
    * `--bin-count` (int): Number of bins when splitting input molecules (default: `50`).
    * `--join-tile-size` (float): Tile size (µm) for molecule–pixel joining (default: `500`).
    * `--max-tile-bytes` (int): Maximum allowed tile size in bytes for PMTiles (default: `5_000_000`).
    * `--max-feature-counts` (int): Maximum number of features per tile (default: `500_000`).
    * `--preserve-point-density-thres` (int): Threshold to preserve point density in PMTiles (default: `1024`).
    * `--transparent-below` / `--transparent-above` (int): Make pixels below/above the threshold transparent for dark/light backgrounds.
    * `--sge-scale` (float): Scales input coordinates to pixels in the output image (default: 1).
    * `--umap-colname-factor` (str): Column name encoding the dominant factor assignment in a UMAP TSV (default: `topK`).
    * `--umap-colname-x` (str): Column name for the UMAP X coordinate (default: `UMAP1`).
    * `--umap-colname-y` (str): Column name for the UMAP Y coordinate (default: `UMAP2`).
    * `--umap-min-zoom` (int): Minimum zoom for generated UMAP PMTiles (default: `0`).
    * `--umap-max-zoom` (int): Maximum zoom for generated UMAP PMTiles (default: `18`).
    * `--keep-intermediate-files` (flag): Retain intermediate files generated.
    * `--skip-raster` (flag): Skip raster tile generation and related dependencies.
    * `--tmp-dir` (str): Path to a temporary directory (default: `<out-dir>/tmp`).

    **Environment Parameters**:

    * `--gzip` (str): Path to the `gzip` binary. For faster compression, use `pigz -p4` (default: `gzip`).
    * `--pmtiles` (str): Path to the `pmtiles` binary from go-pmtiles (default: `pmtiles`).
    * `--gdal_translate` (str): Path to the `gdal_translate` binary (default: `gdal_translate`).
    * `--gdaladdo` (str): Path to the `gdaladdo` binary (default: `gdaladdo`).
    * `--tippecanoe` (str): Path to the `tippecanoe` binary (default: `tippecanoe`).
    * `--spatula` (str): Path to the `spatula` binary (default: `spatula`).

    **Run Parameters**:

    * `--dry-run` (flag): Generate the Makefile; do not execute.
    * `--restart` (flag): Ignore existing outputs and rerun all steps.
    * `--makefn` (str): Makefile name to write (default: `run_cartload2.mk`).
    * `--n-jobs` (int): Number of parallel jobs (default: 1).
    * `--threads` (int): Max threads per job (tippecanoe, etc.) (default: 4).
    * `--log` (flag): Write logs to a file under the output directory.
    * `--log-suffix` (str): Log filename suffix (default: `.log`).

---
## Output

### Copied FICTURE Output

Copied FICTURE output from `<fic_dir>`. See formats in [FICTURE analysis](./run_ficture2.md).

### Rasterized Transcript-level SGE
* SGE mono PMTiles (`sge-mono-dark.pmtiles` and `sge-mono-light.pmtiles`): Rasterized gene expression tiles created from raw SGE for web visualization.

### Rasterized Spatial Factor Maps
* Factor probability PMTiles (`*-results.pmtiles`): Vector tiles encoding posterior topic probabilities per spatial location.
* Decoded factor PMTiles (`*-pixel-raster.pmtiles`): Rasterized spatial factor maps derived from pixel-level decoding results.

### Joined Molecule-factor PMTiles
* Joined molecule-factor TSV (`transcripts_pixel_joined.tsv.gz`): Merged file linking transcript-level SGE with decoded pixel factors.
* Final molecule PMTiles (`genes_bin*.pmtiles`, `genes_index.tsv`, `genes_pmtiles_index.tsv`, `genes_bin_counts.json`): Indexed, multi-feature PMTiles built from joined pixel-factor data for CartoScope.

### Summary Files
* FICTURE assets file (`ficture_assets.json`): JSON catalog listing all output files and their roles for each trained model.
* Catalog file (`catalog.yaml`): Final YAML file summarizing all visual assets and layers for further deployment and visualization.
