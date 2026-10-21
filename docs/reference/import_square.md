# Square Analysis Import

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) calls `import_visiumhd_square` for the Visium HD square-bin analyses (e.g., 8 µm and 16 µm) found under each sample's Space Ranger output. This page is for importing square bins directly, step by step.

## Overview

`import_visiumhd_square` packages one bin size of Space Ranger square-bin analysis outputs (positions, clusters, DE, and UMAP) into a square-polygon PMTiles layer and a square assets JSON for [`run_cartload2 --square-assets`](./run_cartload2.md), and can append the layer to an existing `catalog.yaml`.

---
## Requirements

- Space Ranger square-bin outputs for one bin size (8 or 16 µm are typical), all **required**:
    - Positions Parquet (`spatial/tissue_positions.parquet`)
    - Cluster CSV (`analysis/clustering/gene_expression_graphclust/clusters.csv`)
    - UMAP projection CSV (`analysis/umap/gene_expression_2_components/projection.csv`)
    - Bin-feature MEX (`filtered_feature_bc_matrix/`) to regenerate pseudobulk and DE (default), or the DE CSV (`analysis/diffexp/gene_expression_graphclust/differential_expression.csv`) with `--skip-mex-diffexp`
    - Scale factors JSON (`spatial/scalefactors_json.json`), unless you pass `--scale-json ""` and `--units-per-um`
- Tools (defaults point to the bundled submodules): `tippecanoe`, `spatula`, punkst (`--ficture2`, for the factor report of regenerated DE); `pigz` (or `parquet-tools` with `--use-parquet-tools`); `Rscript` on `PATH` for UMAP plots.

---
## Actions

Converts one bin size into PMTiles (each bin drawn as a `--bin-size` µm square centered on its position), regenerates the cluster pseudobulk and DE from the bin MEX (default), converts the UMAP projection, writes `<outprefix>_assets.json`, and, with `--update-catalog`, appends the layer to an existing `catalog.yaml`.

- `--skip-mex-diffexp`: read DE from the DE CSV instead of regenerating it.
- `--skip-mex-pseudobulk`: skip the pseudobulk; must be combined with `--skip-mex-diffexp` (then no MEX is needed).

---
## Example Usage

=== "JSON mode (auto paths)"
    ```bash
    cartloader import_visiumhd_square \
      --in-json /path/to/space_ranger_assets.json \
      --bin-size 8 \
      --outprefix /path/to/cartload2/spaceranger-sq008 \
      --threads 8
    ```

    Reads the `GRID_<bin-size>um` section of a [`load_space_ranger`](./visiumhd_addon.md) manifest (`POSITION`, `CLUSTER`, `DE`, `UMAP_PROJ`, `GRID_FEATURE_MEX`, `SCALE`).

=== "Manual mode (explicit files)"
    ```bash
    cartloader import_visiumhd_square \
      --in-dir /path/to/space_ranger_output/binned_outputs/square_008um \
      --bin-size 8 \
      --outprefix /path/to/cartload2/spaceranger-sq008 \
      --threads 8
    ```

    File locations default to the layout of one `square_<NNN>um` directory; override `--parquet`, `--csv-clust`, `--csv-diffexp`, `--csv-umap`, `--mex-dir`, or `--scale-json` as needed.

---
## Parameters

### Input/Output
- `--in-json` (str): Space Ranger assets JSON (from `load_space_ranger`); mutually exclusive with `--in-dir`.
- `--outprefix` (str, required): Output prefix for generated files.
- `--bin-size` (int, required): Bin size in µm (e.g., 8 or 16).
- `--id`, `--name` (str): Asset ID and display name (default: basename of `--outprefix`).

### Manual Inputs (relative to `--in-dir`)
- `--in-dir` (str): Square-bin output directory for manual mode.
- `--parquet` (str): Positions Parquet (default: `spatial/tissue_positions.parquet`).
- `--csv-clust` (str): Cluster CSV (default: `analysis/clustering/gene_expression_graphclust/clusters.csv`).
- `--csv-diffexp` (str): DE CSV, used with `--skip-mex-diffexp` (default: `analysis/diffexp/gene_expression_graphclust/differential_expression.csv`).
- `--csv-umap` (str): UMAP projection CSV (default: `analysis/umap/gene_expression_2_components/projection.csv`).
- `--mex-dir` (str): Bin-feature MEX directory (default: `filtered_feature_bc_matrix`); file names via `--mex-bcd`, `--mex-ftr`, `--mex-mtx` (default: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`).
- `--scale-json` (str): Scale factors JSON (default: `spatial/scalefactors_json.json`). It must exist and overrides `--units-per-um`; pass `--scale-json ""` to use `--units-per-um`.
- `--units-per-um` (float): Coordinate units per µm, used only without a scale JSON (default: `1`).

### Catalog
- `--update-catalog` (flag): Append the square layer to `assets.factors` of an existing catalog (the catalog must already contain `assets.factors`; re-running appends a second entry).
- `--catalog-yaml` (str): Catalog to update (default: `catalog.yaml` in the directory of `--outprefix`).
- `--out-catalog-yaml` (str): Write the updated catalog here instead of overwriting `--catalog-yaml`.

??? note "Auxiliary Parameters"
    Recommend to use the default values; override only if needed.

    **PMTiles Conversion**:

    * `--min-zoom` / `--max-zoom` (int): Zoom range for the square layer (default: `10` / `18`).
    * `--umap-min-zoom` / `--umap-max-zoom` (int): Zoom range for UMAP (default: `0` / `18`).
    * `--max-tile-bytes` (int): Max bytes per tile (default: `5000000`).
    * `--max-feature-counts` (int): Max features per tile (default: `500000`).
    * `--preserve-point-density-thres` (int): Density preservation threshold (default: `1024`).

    **Column Names**:

    * `--pos-colname-barcode` / `--pos-colname-x` / `--pos-colname-y` (str): Positions Parquet columns (default: `barcode` / `pxl_col_in_fullres` / `pxl_row_in_fullres`).
    * `--clust-colname-barcode` / `--clust-colname-cluster` (str): Cluster CSV columns (default: `Barcode` / `Cluster`).
    * `--umap-colname-barcode` / `--umap-colname-x` / `--umap-colname-y` (str): UMAP CSV columns (default: `Barcode` / `UMAP-1` / `UMAP-2`).

    **Coloring and DE**:

    * `--tsv-cmap` (str): Cluster color table TSV (default: `<cartloader_dir>/assets/default_color_map.tsv`).
    * `--de-max-pval` (float): Max p-value (default: `0.01`).
    * `--de-min-fc` (float): Min fold change (default: `1.2`).
    * `--skip-mex-diffexp` / `--skip-mex-pseudobulk` (flag): See [Actions](#actions).

    **Run Parameters**:

    * `--threads` (int): Threads for tippecanoe (default: `4`).
    * `--tmp-dir` (str): Temporary directory (default: `tmp/` in the directory of `--outprefix`).
    * `--keep-intermediate-files` (flag): Keep intermediate files.
    * `--log` (flag) / `--log-suffix` (str): Write a log next to the outputs (default suffix: `.log`).

    **Environment Parameters**:

    * `--use-parquet-tools` (flag): Convert the Parquet with `parquet-tools csv` instead of `cartloader parquet_to_csv_rapid` (polars + pigz).
    * `--parquet-tools` (str): `parquet-tools` binary (default: `parquet-tools`).
    * `--pigz` (str) / `--pigz-threads` (int): pigz binary and threads (default: `pigz` / `4`).
    * `--tippecanoe` (str): Default `<cartloader_dir>/submodules/tippecanoe/tippecanoe`.
    * `--spatula` (str): Default `<cartloader_dir>/submodules/spatula/bin/spatula`.
    * `--ficture2` (str): punkst repository (default: `<cartloader_dir>/submodules/punkst`).
    * `--python`, `--gzip`, `--sort` (str): Binaries on `PATH` (default: `python3`, `gzip`, `sort`).
    * `--R` (str): Checked before UMAP plotting (default: `R`); the plots themselves call `Rscript`.

---
## Output

Written next to `--outprefix`:

- `<outprefix>.pmtiles`: Square bins with `barcode` and `topK` (cluster index, ranked by cluster size).
- `<outprefix>-rgb.tsv`: Cluster color table.
- `<outprefix>-bulk-de.tsv`: Cluster DE table.
- `<outprefix>-pseudobulk.tsv.gz` and `<outprefix>-info.tsv`: Cluster pseudobulk and factor info (when regenerated).
- `<outprefix>-umap.tsv.gz`, `<outprefix>-umap.pmtiles`, `<outprefix>.umap.png`, `<outprefix>.umap.single.binary.png`: UMAP outputs.
- `<outprefix>_assets.json`: Square asset entry for `run_cartload2 --square-assets`.

Intermediate files (`.parquet.csv.gz`, `.geojson`, `-clust.tsv.gz`, `.sptsv.*`) are removed unless `--keep-intermediate-files`.
