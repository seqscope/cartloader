# Multi-Sample Cell Analysis (`run_ficture2_multi_cells`)

## Overview

`run_ficture2_multi_cells` runs a cell-level analysis on top of a [`run_ficture2_multi`](./run_ficture2_multi.md) output directory. It builds per-cell count matrices from cell-assigned transcripts or cell-by-gene MEX files, factorizes the cells with LDA (de novo, or by projecting an existing model), clusters them with Leiden, computes UMAP (and optionally t-SNE) manifolds, cluster pseudobulk, differential expression and a factor-by-cluster heatmap, and pixel-decodes each sample with the cluster pseudobulk as the model. Every step runs both jointly across samples and per sample, and the results are recorded in manifests that [`run_cartload2`](./run_cartload2.md) and [`run_cartload2_multi`](./run_cartload2_multi.md) package as cell layers.

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) calls `run_ficture2_multi_cells` for each cell analysis (with `--all`, `--out-prefix <analysis id>`, `--pretrained-model` set to the pixel FICTURE model, and generated `--list-*` files; see [`cell_analyses`](./run_together_inputs.md#cell_analyses)). Use this page to run it directly in step-by-step or custom workflows.

---
## Requirements

- A `run_ficture2_multi` output directory (`--in-dir`, default `--out-dir`) with `samples/<sample>/<sample>.tiled.tsv` and `.tiled.index` for every sample (a `--prepare-only` run is enough).
- A cell-count source for each sample, chosen per sample in this order:
    1. `--mex-list`: a cell-by-gene MEX (`spatula mex2sptsv`);
    2. `--tsv-list`: an external pixel TSV with a header naming `X`, `Y`, `gene`, `count`, `cell_id` (e.g. a Stereo-seq cell-bin file);
    3. otherwise the sample's tiled transcript, whose cell-ID column is `--colidx-cell-id` (default: column 5, as written by `sge_convert --csv-colnames-others cell_id` or `reformat_cosmx`).
- `punkst` (`--ficture2`, default `<cartloader_dir>/submodules/punkst`, with `ext/py/factor_report.py`), `spatula` (default `<cartloader_dir>/submodules/spatula/bin/spatula`), `Rscript` (`--R`), `python3`, `gzip`, and `sort`.
- `cartloader` on the `PATH` (the Makefile calls `cartloader lda_leiden_cluster_fast` and `cartloader lda_tsne`).

---
## Actions

!!! warning "Action Specifications"
    No action runs by default. `--all` enables `--sptsv`, `--lda`, `--leiden`, `--umap`, `--pseudobulk`, `--heatmap`, and `--decode`; `--tsne` must be added explicitly.

- **`--sptsv`**: Builds each sample's cell-by-feature matrix (SPTSV) from its cell-count source, applying the [feature filters](#feature-filters), `--min-feature-count` and (for transcript sources) `--min-cell-count`; then merges the samples into `<out-prefix>.sptsv` and randomizes the order.
- **`--lda`**: Trains an LDA model with `--n-factor` factors on the merged matrix, or projects `--pretrained-model` (plain or `.gz`) onto it; then projects the model onto each sample.
- **`--leiden`**: Leiden clustering of the joint LDA results (`--leiden-resolution`), split per sample; with `--list-cluster`, existing cluster labels are imported instead. Also draws each sample's clusters in space when cell coordinates are available (see [Cell coordinates](#cell-coordinates)).
- **`--umap`** / **`--tsne`**: Shared and per-sample manifolds of the LDA results, colored by cluster.
- **`--pseudobulk`**: Cluster-by-feature pseudobulk matrices (`spatula sptsv2model`), shared and per sample, with a cluster color map.
- **`--heatmap`**: Differential expression of the LDA model and of the cluster pseudobulk, a cluster factor report, and a factor-by-cluster heatmap, shared and per sample. The cluster DE and factor-info files recorded in the manifests are written by this step.
- **`--decode`**: Pixel-level decoding of each sample's tiled transcripts using the shared cluster pseudobulk as the model, with a PNG image.

Step dependencies: `--leiden` needs the LDA; `--umap`/`--tsne` need the LDA and clusters; `--pseudobulk` needs the SPTSV and clusters; `--heatmap` needs the LDA, clusters and pseudobulk; `--decode` needs the pseudobulk. Steps completed in an earlier invocation are recognized by their `.done` files. A run without `--sptsv` must pass `--sptsv-prefix`.

!!! warning "Run the steps together"
    The manifests are rewritten on every invocation and list only the outputs of the steps enabled in that invocation. Use `--all` (plus `--tsne` if wanted), or repeat the earlier step flags when adding a step; completed steps are skipped.

---
## Example Usage

Project an existing pixel model onto cells whose IDs are in the tiled transcripts, with cell coordinates and boundaries from list files (what `run_together` does):

```bash
cartloader run_ficture2_multi_cells \
    --all \
    --out-dir /path/to/out/ficture_multi \
    --out-prefix cells \
    --pretrained-model /path/to/out/ficture_multi/t12_f24.model.tsv \
    --list-xy /path/to/list_xy.tsv \
    --list-boundaries /path/to/list_boundaries.tsv \
    --threads 8 \
    --n-jobs 4
```

Train a de-novo cell model from cell-by-gene MEX files:

```bash
cartloader run_ficture2_multi_cells \
    --all \
    --out-dir /path/to/out/ficture_multi \
    --out-prefix mexcells \
    --n-factor 24 \
    --mex-list /path/to/list_mex.tsv
```

---
## Input lists

All list files are tab-delimited, one sample per line.

| Flag | Line format | Notes |
|------|-------------|-------|
| `--list-samples` | `SAMPLE_ID ...` | Only the first column is used; each sample needs `<in-dir>/samples/<SAMPLE_ID>/`. Default: every directory under `<in-dir>/samples/`. |
| `--mex-list` | `SAMPLE_ID MEX_DIR`, or `SAMPLE_ID BCD FTR MTX` | With a directory, file names come from `--mex-bcd`/`--mex-ftr`/`--mex-mtx`. |
| `--tsv-list` | `SAMPLE_ID PIXEL_TSV` | The TSV has a header with `X`, `Y`, `gene`, `count`, `cell_id`. A sample in both `--mex-list` and `--tsv-list` uses the MEX. |
| `--list-cluster` | `SAMPLE_ID CLUSTER_FILE` | Each cluster file has two columns, cell ID and cluster ID (tab, comma or space; the first line is skipped as a header unless its cluster ID is a number). IDs are 1-based unless `--zero-based-clust-id`. Samples not listed get a placeholder file. |
| `--list-xy` | `SAMPLE_ID XY_FILE` | Cell coordinates in µm; columns named by `--xy-colname-cell-id`/`-x`/`-y` (tab, comma or space). |
| `--list-boundaries` | `SAMPLE_ID BOUNDARIES_FILE [SCALE_JSON]` | Cell boundary polygons. Recorded in the manifest for packaging, or, with `--boundaries-format geojson`, used to derive cell centroids. |

### Cell coordinates

The per-sample cluster scatter (`.leiden.xy.*`, the source of the cell-point layer) uses, in order: the sample's `--list-xy` file; centroids derived from its `--list-boundaries` GeoJSON when `--boundaries-format geojson` is set (coordinates divided by the `--boundaries-units-key` value of the optional `SCALE_JSON`, cell IDs formatted with `--boundaries-cell-id-format`); or the cell metadata written by `pixel2sptsv` for transcript-based samples. A MEX sample with none of these has no coordinates, and its scatter is skipped.

---
## Feature filters

`--include-feature-list`, `--exclude-feature-list`, `--include-feature-regex`, and `--exclude-feature-regex` are applied when the per-cell matrices are built (`--sptsv`), so a dropped feature is absent from the whole cell analysis — LDA, clusters, pseudobulk, DE, and decode. The tiled transcripts are not changed.

- Regexes use **substring** semantics (`spatula`, `std::regex_search`).
- `spatula` takes at most one include-type and one exclude-type filter, so a list and a regex of the same polarity are first folded into one list with [`feature_select`](./feature_customization.md#resolving-filters-into-one-list-feature_select): an include pair becomes `<out-dir>/<out-prefix>.include_features.tsv` (the list filtered by the regex); an exclude pair becomes `<out-dir>/<out-prefix>.exclude_features.tsv` (the list plus the features matching the regex in `<in-dir>/multi.union_features.tsv`). These files are written when the Makefile is generated.
- A run that reuses `--sptsv-prefix` inherits the filters of the run that built it.

`run_together` passes only its **ingest** filters here; the FICTURE analysis restriction reaches the cell LDA through the projected pixel model instead, so the cell pseudobulk and DE keep every gene (see [Feature filtering: two independent layers](./run_together_inputs.md#feature-filtering-two-independent-layers)).

---
## Parameters

### Action Flags
- `--all`, `--sptsv`, `--lda`, `--leiden`, `--umap`, `--tsne`, `--pseudobulk`, `--heatmap`, `--decode`: See [Actions](#actions).

### Input/Output
- `--out-dir` (str, required): Output directory.
- `--in-dir` (str): The `run_ficture2_multi` output directory (default: `--out-dir`).
- `--out-prefix` (str): Prefix of all outputs and the cell analysis ID (default: `cells`).
- `--list-samples`, `--mex-list`, `--tsv-list`, `--list-cluster`, `--list-xy`, `--list-boundaries` (str): See [Input lists](#input-lists).
- `--mex-bcd`, `--mex-ftr`, `--mex-mtx` (str): MEX file names inside a `--mex-list` directory (defaults: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`).
- `--sptsv-prefix` (str): Prefix of an existing merged SPTSV, required when `--sptsv` is not set and not allowed with it. Normally `<out-dir>/<out-prefix>.sptsv` from an earlier `--sptsv` run (the per-sample matrices and `--pseudobulk` always read from there).
- `--out-json` (str): Accepted but currently unused; the manifests are written as described in [Output](#output).

### Key Parameters
- `--n-factor` (int): Number of LDA factors; required with `--lda` unless `--pretrained-model` is given.
- `--pretrained-model` (str): Model TSV (plain or `.gz`) to project instead of training.
- `--leiden-resolution` (float): Leiden resolution (default: 1.0).
- `--cmap-file` (str): Color map TSV (default: `<cartloader_dir>/assets/default_color_map.tsv`).
- `--decode-fit-width` (int): Hexagon width (µm) for decoding (default: 18).
- `--anchor-resolution` (int): Anchor resolution; the decode uses `fit_width // anchor_resolution + 1` sliding moves (default: 6).

??? note "Auxiliary Parameters"

    **Cell matrix**:

    - `--colidx-x`, `--colidx-y`, `--colidx-feature`, `--colidx-count`, `--colidx-cell-id` (int, 1-based): Columns of the tiled transcript (defaults: 1, 2, 3, 4, 5).
    - `--ignore-ids` (str): Comma-separated cell IDs treated as unassigned (default: `UNASSIGNED,NA,0,-1`).
    - `--min-feature-count` (int): Minimum count per feature for the cell matrices and the pseudobulk (default: 20).
    - `--min-cell-count` (int): Minimum count per cell (transcript sources only; default: 50).

    **Clusters and coordinates**:

    - `--zero-based-clust-id` (flag): Cluster IDs in `--list-cluster` files are already 0-based.
    - `--xy-colname-cell-id`, `--xy-colname-x`, `--xy-colname-y` (str): Columns of the `--list-xy` files (defaults: `cell_id`, `X`, `Y`).
    - `--boundaries-format` (str): `geojson` to derive centroids from `--list-boundaries` files (default: none, boundaries are pass-through).
    - `--boundaries-cell-id-format` (str): Python format string applied to each GeoJSON cell ID so it matches the clustering barcodes (default: `cellid_{:09d}-1`, the Visium HD convention; `{}` keeps the raw ID).
    - `--boundaries-cell-id-prop` (str): GeoJSON property holding the cell ID (default: `cell_id`).
    - `--boundaries-units-key` (str): Key in the `SCALE_JSON` giving microns per pixel (default: `microns_per_pixel`).

    **LDA, DE and decoding**:

    - `--seed` (int): Random seed (default: 1).
    - `--train-epoch` (int): Accepted but currently has no effect (the LDA always runs 2 epochs).
    - `--de-min-ct-per-feature` (int): Minimum count per feature for DE (default: 20).
    - `--de-max-pval` (float): P-value cutoff for DE (default: 1e-3).
    - `--de-min-fold` (float): Fold-change cutoff for DE (default: 1.5).
    - `--decode-pixel-res` (float): Pixel resolution (µm) of the decode (default: 0.5).
    - `--single-molecule` (flag): Decode each molecule instead of pixels.
    - `--decode-scale` (int): Scale from coordinates to pixels of the decode PNG (default: 1).

    **Feature filters**: `--include-feature-list`, `--exclude-feature-list`, `--include-feature-regex`, `--exclude-feature-regex` (see [Feature filters](#feature-filters)).

    **Environment**:

    - `--ficture2` (str): Path to the `punkst` repository (default: `<cartloader_dir>/submodules/punkst`).
    - `--spatula` (str): Path to `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`).
    - `--R` (str): Command that runs the R scripts (default: `Rscript`).
    - `--python` (str): Python 3 for `punkst`'s factor report (default: `python3`).
    - `--gzip` (str), `--sort` (str): Paths to `gzip` and `sort` (defaults: `gzip`, `sort`).

    **Run**:

    - `--dry-run` (flag): Write the Makefile and print the commands only.
    - `--restart` (flag): Ignore existing outputs and rerun all steps.
    - `--threads` (int): Max threads per job (default: 8).
    - `--n-jobs` (int): Parallel Makefile jobs (default: 2).
    - `--makefn` (str): Makefile name under `--out-dir` (default: `run_ficture2_multi_cells.<out-prefix>.mk`).

---
## Output

With `<p>` = `--out-prefix`. Shared outputs under `--out-dir`:

- `<p>.sptsv.tsv`, `.json`, `.feature.counts.tsv`, `.randomized.tsv`, `.list.tsv`: Merged cell matrix.
- `<p>.lda.model.tsv`, `<p>.lda.results.tsv`: The cell LDA model and the cells' factor posteriors.
- `<p>.leiden.tsv.gz`: Cluster of every cell (`sample_id`, `cell_id`, `topK`).
- `<p>.umap.tsv.gz`, `<p>.umap.png`, `<p>.umap.leiden.tsv.gz` (and `<p>.tsne.*` with `--tsne`): Shared manifolds.
- `<p>.leiden.pseudobulk.tsv`, `.cmap.tsv`, `.de.tsv`, `.factor.info.tsv`, `.factor.info.html`: Cluster pseudobulk, colors, DE and report.
- `<p>.lda.model.de.tsv`, `<p>.heatmap.pdf`, `<p>.heatmap.normfrac.tsv`: Model DE and the factor-by-cluster heatmap.
- `ficture.multi.<p>.params.json`: Shared manifest — the per-sample manifests (relative paths) and the shared model, cluster pseudobulk/DE/info, heatmap and manifolds of the steps that ran.
- `.done` flags for each step.

Per sample, in `samples/<sample>/` with prefix `<sample>.<p>`:

- `.sptsv.tsv`, `.sptsv.json`, `.sptsv.feature.counts.tsv`: The sample's cell matrix (plus `.sptsv.cell.metadata.tsv` or `.sptsv.cell.xy.tsv` with cell coordinates).
- `.lda.results.tsv`: The sample's factor posteriors.
- `.leiden.tsv.gz`: Cluster of every cell; `.leiden.xy.tsv.gz` and `.leiden.xy.png` when coordinates are available.
- `.umap.*` (and `.tsne.*`): Per-sample manifolds.
- `.leiden.pseudobulk.tsv`, `.cmap.tsv`, `.de.tsv`, `.factor.info.tsv`: Per-sample cluster pseudobulk, colors, DE and report.
- `.heatmap.pdf`, `.heatmap.normfrac.tsv`: Per-sample heatmap.
- `.pixel.bin` (with its index) and `.pixel.png`: Pixel-level decode.
- `ficture.<p>.params.json`: Per-sample manifest, a `cell_params` object passed to `run_cartload2 --in-cell-params` (auto-detected by `run_cartload2_multi`). Its keys follow the steps that ran: `model_path`/`fit_path` (`--lda`), `sptsv_prefix` (`--sptsv`), `cluster_path` and `cell_xy_path` (`--leiden`), `cmap`/`cluster_pseudobulk`/`cluster_de`/`cluster_info` and their `shared_*` versions (`--pseudobulk`), `cluster_model_heatmap_*` (`--heatmap`), `pixel_bin_prefix`/`pixel_png_path`/`pixel_res`/`decode_scale` (`--decode`), `manifolds`, `analysis_type` (multi-sample runs), and `cell_boundaries_path` (pass-through boundaries).

Both manifests are written when the Makefile is generated.

??? note "Temporary files"
    `lda4hex` uses `<out-dir>/<p>.lda_lda4hex` and `samples/<sample>/<sample>.<p>.lda_lda4hex`; `pixel-decode` uses `samples/<sample>/<sample>.<p>.pixel.tmp`. The R scripts use the system temporary directory.
