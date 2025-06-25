# Spatial Asset Packaging

## Overview

Following spatial factor inference via FICTURE analysis, `cartloader` offers the `run_cartload2` module package SGE data and spatial factor output from [FICTURE analysis](./run_ficture2.md) into standardized, spatially indexed, and storage-efficient [PMTiles](https://github.com/protomaps/PMTiles), a web-native tiling format. These PMTiles outputs are optimized for downstream analysis, interactive web visualization (e.g., in CartoScope), and data sharing across platforms.

## Requirements

- A completed FICTURE run (pixel-level decoding outputs) (from `run_ficture2`)
- A metadata file describing the input-output structure (`ficture.params.json`)
- Pre-installed CLI tools: `tippecanoe`, `gdal_translate`, `gdaladdo`, `pmtiles`, `spatula`, `gzip`

## Example Usage

```bash
cartloader run_cartload2 \
    --makefn run_cartload2.mk \
    --fic-dir /path/to/run_ficture2/results \
    --out-dir /path/to/output/directory \
    --id dataset_id \                            ## replace dataset_id with the id for this dataset
    --colname-count count \
    --n-jobs 20  \
    --threads 20 \
    --spatula /path/to/spatula/binary \
    --pmtiles /path/to/pmtiles/binary \
    --tippecanoe /path/to/tippecanoe/binary
```

## Actions

Specifically, it will perform all the following steps:

- Converts transcript-level SGE to raster-format PMTiles.
- Converts topic proportions (.results.tsv.gz) into vector PMTiles for spatial factors.
- Processes each decoded pixel-level output to generate raster PMTiles for visual overlays.
- Create a joined molecule-feature matrix by joining decoded pixel-level spatial factors from FICTURE with transcript-level molecules from the original SGE based on spatial proximity.
- Converts the joined molecule-feature matrix into storage-efficient, multi-feature PMTiles.
- Generates a JSON file listing all FICTURE assets and a YAML catalog describing the final visualization layers.

## Parameters

The following outlines the **minimum required parameters** for running spatial asset packaging. 

For auxiliary parameters, we recommend using the default values unless you possess a thorough understanding of `run_cartload2`. For further details, refer to the collapsible sections below or run:

```bash
cartloader run_cartload2 --help
```

### Input/Output Parameters

* `--fic-dir` (str): Path to the input directory containing FICTURE output.
* `--out-dir` (str): Path to the output directory for storing generated assets.

### Dataset ID and Descriptions

* `--id` (str): Unique identifier for the output asset set.
* `--title` (str): Optional title for the output assets.
* `--desc` (str): Optional description for the output assets.

??? note "Auxiliary `run_cartload2` Paramaters"

    **Auxiliary Conversion Parameters**:
    
    * `--in-fic-params` (str): Path to input JSON file with SGE files and FICTURE parameters (Default: `ficture.params.json`).
    * `--out-fic-assets` (str): Path to output JSON file to write FICTURE assets (Default: `ficture_assets.json`).
    * `--out-catalog` (str): Path to output YAML file for assets (Default: catalog.yaml).
    * `--background-assets` (str list): List of background asset descriptors in `[id:file]` or `[id1:id2:file]` format.
    * `--rename-x` (str): Column renaming rule for X axis in `tippecanoe` (Default: x:lon).
    * `--rename-y` (str): Column renaming rule for Y axis in `tippecanoe` (Default: y:lat).
    * `--colname-feature` (str): Column name for gene/feature name (Default: gene).
    * `--colname-count` (str): Column name for feature count (Default: count).
    * `--out-molecules-id` (str): Prefix for output molecule PMTiles files (Default: genes).
    * `--max-join-dist-um` (float): Maximum join distance (µm) between molecules and pixels (Default: 0.1).
    * `--join-tile-size` (float): Tile size (µm) for molecule–pixel joining. (Default: 500).
    * `--max-tile-bytes` (int): Maximum allowed tile size in bytes for PMTiles (Default: 5e6).
    * `--max-feature-counts` (int): Maximum number of features per tile (Default: 5e5).
    * `--preserve-point-density-thres` (int): Threshold to preserve point density in PMTiles (Default: 1024).
    * `--keep-intermediate-files` (flag): If set, retain intermediate files generated.
    * `--skip-raster` (flag): If set, skip raster tile generation and related dependencies.
    * `--tmp-dir` (str): Path to a temporary directory (Default: `<out-dir>/tmp`).

    **Auxiliary Environment Parameters**:

    * `--gzip` (str): Path to the `gzip` binary. For faster compression, use `"pigz -p4"` (Default: gzip).
    * `--pmtiles` (str): Path to the `pmtiles` binary from go-pmtiles (Default: pmtiles).
    * `--gdal_translate` (str): Path to the `gdal_translate` binary (Default: gdal_translate).
    * `--gdaladdo` (str): Path to the `gdaladdo` binary (Default: gdaladdo).
    * `--tippecanoe` (str): Path to the `tippecanoe` binary (Default: tippecanoe).
    * `--spatula` (str): Path to the `spatula` binary (Default: spatula).

## Output

### Copied FICTURE Output

Copied FICTURE output from `<fic_dir>`. See formats in [FICTURE analysis](./run_ficture2.md).

### Rasterized transcript-level SGE
* SGE mono PMTiles (`sge-mono-dark.pmtiles` and `sge-mono-light.pmtiles`): Rasterized gene expression tiles created from raw SGE for web visualization.

### Rasterized Spatial Factor Maps
* Factor probability PMTiles (`*-results.pmtiles`): Vector tiles encoding posterior topic probabilities per spatial location.
* Decoded factor PMTiles (`*-pixel-raster.pmtiles`): Rasterized spatial factor maps derived from pixel-level decoding results.

### Joined molecule-factor PMTiles
* Joined molecule-factor TSV (`transcripts_pixel_joined.tsv.gz`): Merged file linking transcript-level SGE with decoded pixel factors.
* Final molecule PMTiles (*_pmtiles_index.tsv, `*_bin_counts.json`): Indexed, multi-feature PMTiles built from joined pixel-factor data for CartoScope.

### Summary Files
* FICTURE assets file (ficture_assets.json): JSON catalog listing all output files and their roles for each trained model.
* Catalog file (catalog.yaml): Final YAML file summarizing all visual assets and layers for further deployment and visualization.