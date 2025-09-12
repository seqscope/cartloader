# Multi‑Sample Spatial Asset Packaging

## Overview

`run_cartload2_multi` batches CartLoad2 packaging across multiple samples. This module is specifically designed for packaging jointly FICTURE runs ([`run_ficture2_multi`](./run_ficture2_multi.md)), which generates per‑sample FICTURE assets and want per‑sample PMTiles and catalogs.

For each sample, it invokes `run_cartload2` to package SGE and optional FICTURE outputs into PMTiles and a catalog. It generates a Makefile to parallelize work and ensures all per‑sample catalogs are produced consistently.

!!! warning "Scope and Limitations"
    Convenience wrapper to package SGE and FICTURE outputs produced by `run_ficture2_multi`.
    It does not handle cell analysis results or background images (e.g., histology).
    To include all assets, run `run_cartload2` per sample or pair `run_ficture2_multi` with `import_image` and `import_*_cell`.

---
## Requirements
- TSV file containing the information for input samples (two to four columns: `id`, `transcript_path`, optional `title`, optional `desc`)  
- Per‑sample FICTURE outputs under `<fic_dir>/samples/<sample_id>/` with `ficture.params.json` (from `run_ficture2_multi`).
- Unified SGE and any deployed images/cells (as required by each sample’s `run_cartload2`).
- CLI tools used by CartLoad2: `tippecanoe`, `gdal`, `pmtiles`, `gzip`/`pigz`, `spatula`.

---
## Example Usage

```bash
cartloader run_cartload2_multi \
  --in-list /path/to/samples.tsv \
  --fic-dir /path/to/run_ficture2_multi/out \
  --out-dir /path/to/out/cartload2_multi \
  --spatula /path/to/spatula/binary \
  --pmtiles /path/to/pmtiles/binary \
  --tippecanoe /path/to/tippecanoe/binary \
  --n-jobs 10 \
  --threads 10
```

---
## Actions

For each row in `--in-list`, apply:

- Locate per‑sample FICTURE assets at `--fic-dir/samples/<id>/ficture.params.json`.
- Run `cartloader run_cartload2` with `--fic-dir` set to the sample subfolder and `--id` normalized from `<id>` (lowercase, `_` → `-`), and save output to `<out-dir>/<id>/`.
- Write a per‑sample `catalog.yaml` (or `--out-catalog` name) and PMTiles under `<out-dir>/<id>/`.
- Emit a top‑level Makefile to run samples in parallel.

---
## Parameters

Below lists the most common parameters. Additional flags are forwarded to per‑sample `run_cartload2`.

### Input/Output Parameters

- `--in-list` (str, required): Two‑ to four‑column TSV: `<sample_id> <transcript.tsv[.gz]> [title] [desc]`.
- `--fic-dir` (str, required): FICTURE output root containing `samples/<id>/ficture.params.json`.
- `--out-dir` (str, required): Output root for per‑sample catalogs and PMTiles.
- `--out-catalog` (str): Catalog filename per sample (default: `catalog.yaml`).
- `--gdal_translate` (str): Path to `gdal_translate` (forwarded to `run_cartload2`).


??? note "Auxiliary Parameters"
    Recommend to use the default values; override only if needed.

    **Auxiliary Conversion Parameters**:
    
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

Outputs are written under `--out-dir`:

- Per‑sample catalogs: `<out-dir>/<id>/catalog.yaml` (or `--out-catalog` name).
- Per‑sample PMTiles and asset JSONs produced by `run_cartload2` under `<out-dir>/<id>/`.
- Makefile: `<out-dir>/<makefn>` capturing all per‑sample targets.

See `run_cartload2.md` for detailed CartLoad2 output formats (SGE rasters, factor PMTiles, joined molecule PMTiles, summaries).

---
## See Also

- Reference: `run_cartload2.md` — Single‑sample packaging and outputs
- Reference: `run_ficture2_multi.md` — Multi‑sample FICTURE pipeline
