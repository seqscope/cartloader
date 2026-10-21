# Spatial Factor Inference Analysis using FICTURE

## Overview
Following format conversion, `CartLoader` provides the **`run_ficture2`** module to run spatial factor inference on a single sample using [**`FICTURE`**](https://www.nature.com/articles/s41592-024-02415-2) (Si et al., *Nature Methods*, 2024), as implemented in `punkst`. This method infers spatial factors directly at the pixel level with submicron resolution, eliminating the need for segmentation.

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) runs FICTURE as part of the pipeline, through [`run_ficture2_multi`](./run_ficture2_multi.md) (also for a single sample). Use `run_ficture2` directly in step-by-step or custom single-sample workflows.

!!! question
    - [What is `FICTURE` and `punkst`?](../faq/ficture.md)
    - [How can I speed up runs safely?](../faq/speed_up_runs.md)
    - [How do I resume a run, or rerun everything from scratch?](../faq/restart_or_resume.md)

---
## Requirements

- An SGE in the unified format (from [SGE format conversion](./sge_convert.md)): a transcript TSV with a header line, X/Y at `--colidx-x`/`--colidx-y`, and the feature name and count in columns **3 and 4** (fixed positions).
- `punkst`: `--ficture2` must point to a `punkst` repository containing `bin/punkst` and `ext/py/factor_report.py` (default: the repository submodule `<cartloader_dir>/submodules/punkst`).
- `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`), `gzip`, `sort`, and `python3`.
- `Rscript` on the `PATH` for `--umap`.
- `cartloader` on the `PATH` (the Makefile calls `cartloader write_json_for_ficture2`).

## Actions

!!! warning "Action Specifications"
    No action runs by default. Activate at least one using the [action parameters](#action-parameters).

### Tiling Step (`--tile`)
The tiling step takes the standardized SGE (from [SGE format conversion](./sge_convert.md)) as input. It reorganizes input coordinate data into non‑overlapping square tiles in a plain TSV format and generates an index file with tile offsets to enable efficient random access.

### Segmentation Step (`--segment`)
Aggregates the tiled pixel data into non-overlapping hexagons (one width per `--width` value), writing a tab-delimited file of hexagon records and a JSON metadata file. Requires the tiling outputs.

### LDA Training Step (`--init-lda`)
Trains a Latent Dirichlet Allocation (LDA) model on the hexagon counts for each (`--width`, `--n-factor`) pair, then writes the factor color map, a differential-expression table and a factor report.

### Decoding Step (`--decode`)
Applies each trained model to the tiled pixel-level data to infer the top factors and their posterior probabilities at each pixel. It writes binary pixel-level results, a pseudobulk feature-by-factor matrix, a differential-expression table, a factor report, and a PNG image.

### UMAP Visualization Step (`--umap`)
Generates a UMAP embedding of the hexagon-level factor posteriors from [LDA training](#lda-training-step-init-lda), writing both coordinates and plots for reuse downstream.

### 10x Segmentation Step (`--segment-10x`)
Aggregates the tiled pixel data into hexagons in 10x Genomics MEX format (minimum count 1 per hexagon), for use with tools designed for 10x data.

!!! note "Step dependencies"
    `--decode` and `--umap` need `--init-lda` in the same invocation. A model that is already trained is not retrained: completed steps are skipped through their `.done` flags, so rerun with the same flags plus the new action.

---
## Example Usage

```bash
cartloader run_ficture2 \
    --main \
    --in-transcript /path/to/sge/transcripts.unsorted.tsv.gz \
    --in-feature /path/to/sge/feature.clean.tsv.gz \
    --in-minmax /path/to/sge/coordinate_minmax.tsv \
    --out-dir /path/to/output/directory \
    --width 18 \
    --n-factor 24 \
    --exclude-feature-regex '^(mt-.*|Gm[0-9]+)$' \
    --n-jobs 20 \
    --threads 20
```

`--in-json /path/to/sge/sge_assets.json` (written by `sge_convert`) can replace the three `--in-*` flags.

!!! warning "Feature regexes must match the whole name"
    `--include-feature-regex` / `--exclude-feature-regex` are passed to `punkst lda4hex`, which requires the pattern to match the **entire** feature name (`std::regex_match`). `'^(mt-|Gm\d+$)'` would therefore exclude `Gm123` but no `mt-` gene; write `'^(mt-.*|Gm[0-9]+)$'` instead. This differs from `sge_convert`, `run_ficture2_multi` and `feature_filtering`, which match substrings (see [Feature Customization](./feature_customization.md#regex-semantics-differ-between-commands)).

---
## Parameters

Below are the core parameters. See more details in the collapsible section ("Auxiliary Parameters") below.

#### Action Parameters

* `--main`: Run `--tile`, `--segment`, `--init-lda`, `--decode`, and `--umap` (not `--segment-10x`).
* `--tile`: Run [tiling step](#tiling-step-tile).
* `--segment`: Run [segmentation step](#segmentation-step-segment).
* `--init-lda`: Run [LDA training step](#lda-training-step-init-lda).
* `--decode`: Run [decoding step](#decoding-step-decode).
* `--umap`: Run [UMAP step](#umap-visualization-step-umap).
* `--segment-10x`: Run [10x segmentation step](#10x-segmentation-step-segment-10x).

#### Input/Output Parameters

* `--out-dir` (str, required): Output directory to store all result files.
* `--in-transcript` (str): Input transcript-indexed SGE TSV (default: `<out-dir>/transcripts.unsorted.tsv.gz`).
* `--in-minmax` (str): Input coordinate min/max TSV (default: `<out-dir>/transcripts.tiled.coord_range.tsv`, written by tiling). Pass it explicitly so that `ficture.params.json` records it (see [Summarization Output](#summarization-output)).
* `--in-feature` (str): Input per-feature count TSV with a header line; the first column is the feature name and the last the count. Its features define the feature dictionary used by segmentation and LDA. Default: the feature list written by tiling (`transcripts.tiled.features.tsv`).
* `--in-json` (str): Input manifest JSON with `transcript`, `feature`, and `minmax` keys (typically `sge_assets.json` from `sge_convert`); replaces the three flags above.
* `--out-json` (str): Output JSON summarizing FICTURE parameters (default: `<out-dir>/ficture.params.json`).

#### Key Parameters

* `--width` (str): Comma-separated hexagon flat-to-flat widths (µm) for segmentation and LDA training. Required by `--segment`, `--init-lda`, `--decode`, and `--umap`.
* `--n-factor` (str): Comma-separated factor counts for LDA training (default: `12`, with a warning). Every (width, factor count) pair is trained, as model `t{width}_f{n_factor}`.
* `--include-feature-regex` (str): Regex of feature names to include in LDA training (whole-name match).
* `--exclude-feature-regex` (str): Regex of feature names to exclude from LDA training (whole-name match).
* `--cmap-file` (str): Color map TSV; its first `n_factor` colors are used (default: `<cartloader_dir>/assets/default_color_map.tsv`).
* `--segment-width-10x` (str): Comma-separated hexagon widths (µm) in 10x format (required with `--segment-10x`).

??? note "Auxiliary Parameters"

    **Auxiliary Input Parameters**

    * `--in-feature-ficture` (str): A feature file used instead of `--in-feature` for the FICTURE analysis (same format), restricting the analysis to the listed features.
    * `--colidx-x`, `--colidx-y` (int): 1-based column indices of X and Y in the transcript file (defaults: 1, 2).
    * `--colname-feature`, `--colname-count` (str): Header names written to the FICTURE feature file `*.hdr.tsv` (defaults: `gene`, `count`). They do not select input columns.
    * `--model-id` (str): Name of the model in place of `t{width}_f{n_factor}`; only valid when a single width and a single factor count are given.

    **Auxiliary FICTURE Parameters**:

    * Tiling-specific parameters:
        * `--tile-size` (int): Tile size (default: 500).
        * `--tile-buffer` (int): Per-tile buffer size (default: 1000).
    * Segmentation-specific parameters:
        * `--min-ct-per-unit-hexagon` (int): Minimum count per hexagon (default: 50).
    * LDA training-specific parameters:
        * `--min-ct-per-unit-train` (int): Minimum count per hexagon for training (default: 50).
        * `--minibatch-size` (int): Minibatch size (default: 500).
        * `--train-epoch` (int): Number of training epochs (default: 2).
    * Decoding-specific parameters:
        * `--fit-width` (int): Hexagon width (µm) for decoding; give a single value (default: the training width).
        * `--anchor-res` (int): Anchor resolution; the decode uses `fit_width / anchor_res` sliding moves (default: 6).
        * `--decode-pixel-res` (float): Pixel resolution (µm) of the decode (default: 0.5).
        * `--single-molecule` (flag): Decode each molecule instead of pixels at `--decode-pixel-res`.
        * `--decode-scale` (int): Scale from input coordinates to pixels of the output PNG (default: 1).
    * Shared parameters across steps:
        * `--seed` (int): Random seed (default: 1).
        * `--min-ct-per-feature` (int): Minimum count per feature for LDA training and differential expression (default: 20).
        * `--de-max-pval` (float): P-value cutoff for differential expression (default: 1e-3).
        * `--de-min-fold` (float): Fold-change cutoff for differential expression (default: 1.5).

    **Environment Parameters**:

    * `--ficture2` (str): Path to the `punkst` repository (default: `<cartloader_dir>/submodules/punkst`).
    * `--spatula` (str): Path to `spatula` (default: `<cartloader_dir>/submodules/spatula/bin/spatula`).
    * `--gzip` (str): Path to `gzip`; consider `"pigz -p 4"` for speed (default: `gzip`).
    * `--sort` (str): Path to `sort` (default: `sort`).
    * `--sort-mem` (str): Memory per `sort` process (default: `1G`).
    * `--python` (str): Python 3 used for `punkst`'s factor report (default: `python3`).
    * `--R` (str): R binary checked with `--umap` (default: `R`); the UMAP scripts themselves run with `Rscript` from the `PATH`.

    **Run Parameters**:

    * `--dry-run` (flag): Write the Makefile and print the commands without executing.
    * `--restart` (flag): Ignore existing outputs and rerun all steps.
    * `--makefn` (str): Makefile name under `--out-dir` (default: `run_ficture2.mk`).
    * `--n-jobs` (int): Number of parallel jobs (default: 1).
    * `--threads` (int): Max threads per job (default: 4).

??? note "Temporary files"
    `punkst` writes temporary files under `<out-dir>/tmp` (tiling), `<out-dir>/tmp/<width>` (segmentation) and `<out-dir>/tmp/<model_id>_pixel` (decoding). The UMAP R scripts run with `TMPDIR` set to `<out-dir>/<model_id>_create_umap`, `_draw_umap` and `_draw_umap_single` (created and removed per step), so large inputs do not fill the system `/tmp`; the UMAP nearest-neighbor index goes to `<out-dir>/<model_id>.umap.tmp` (removed on completion). `lda4hex` is not given a temporary directory and uses `punkst`'s default.

---
## Output

### Tiling Output

* `transcripts.tiled.tsv`: A tab-delimited TSV file where each line records a pixel with X and Y coordinates, a feature name (e.g., gene), and its associated count, grouped by tile.
    ```text
    3815.69	491.9	Arfgef1 1
    ```
* `transcripts.tiled.index`: An index file that records the byte offsets of each tile in the tiled TSV file, enabling fast random access to tile-specific data. It includes comment lines with metadata (tile size, pixel count, min/max X and Y coordinates) and data lines listing the offset positions for each tile.
    ```text
    # tilesize	500
    # npixels	333021294
    # xmin	0.28
    # xmax	19851.1
    # ymin	0.5
    # ymax	13476.4
    0	7	0	        13103550	569474
    0	8	13103550	40259340	1180529
    ```
* `transcripts.tiled.coord_range.tsv`: A tab-delimited TSV file provides xmin, xmax, ymin, and ymax as key–value pairs.
    ```text
    xmin	0.28
    xmax	19851.1
    ymin	0.5
    ymax	13476.4
    ```
* `transcripts.tiled.features.tsv`: A headerless TSV listing each feature name and its total count.
    ```text
    1810009J06Rik	1
    Gm48545	        1
    ```
* FICTURE feature files: `transcripts.tiled.features.hdr.tsv` (the same list with a header), or, with `--in-feature`/`--in-feature-ficture`, `transcripts.tiled.selected_features.tsv` (headerless) and `transcripts.tiled.selected_features.hdr.tsv`.

### Segmentation Output

* `hexagon.d_{width}.randomized.tsv`: A plain-text TSV file that stores sparse feature counts for hexagon units, one hexagon per line, sorted by a random key.

    !!! info "File Format"
        This is not a regular table—the number of columns varies across lines depending on the number of features present in each hexagon.

    ```text
    00001c64  18999.0000  9773.9627  227  244   2417 1  10469 2  2448 1  543 1  10935 1  6460 1  ...
    ```

    * 1st column (str): Random identifier for the hexagon
    * 2nd and 3rd columns (float): Coordinates of the hexagon center
    * 4th and 5th columns (int): Number of unique features and total counts, repeated per modality (for K modalities)
    * 6th column onward (integer integer): Index–count pairs for non-zero features, using the 0-based indices of the feature dictionary

* `hexagon.d_{width}.json`: The metadata for the hexagons.
    ```json
    {
    "dictionary": {
        "A2m": 9636,
        "AA986860": 6486,
        ...
        },
        "header_info": [
            "random_key",
            "x",
            "y"
        ],
        "hex_size": 10.392304845413264,
        "icol_x": 1,
        "icol_y": 2,
        "n_features": 11319,
        "n_modalities": 1,
        "n_units": 498019,
        "offset_data": 3,
        "random_key": 0
    }
    ```
    * `dictionary`: Feature name–index pairs
    * `header_info`: Names of the leading columns
    * `hex_size`: Hexagon side length (width / √3)
    * `icol_x` and `icol_y`: 0-based column indices of the X and Y coordinates
    * `n_features`: Number of features
    * `n_modalities`: Number of modalities
    * `n_units`: Number of hexagons
    * `offset_data`: Number of leading columns before the per-hexagon counts
    * `random_key`: 0-based column index of the random key
* `hexagon.d_{width}.done`: done flag.

### LDA Training Output

* `t{width}_f{n_factor}.model.tsv`: The learned feature-by-factor matrix.
    ```text
    Feature  0        1      2      3      4      5      6         7      8      9      10        11       12     13     14        15     16     17     18      19     20      21     22     23
    Neu4     481.168  0.042  0.042  0.042  0.042  0.042  1828.293  0.042  0.042  0.042  2290.443  836.722  0.042  0.042  1704.589  0.042  0.042  0.042  29.738  0.042  24.854  0.042  0.042  0.042
    ```
    * Rows (str): Features (e.g., gene names)
    * Columns (int): Factor IDs
    * Values (float): Weight of a feature in each factor
* `t{width}_f{n_factor}.results.tsv.gz`: A tab-delimited file containing the posterior factor distributions for each hexagon.
    ```text
    x           y          0       1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21      22      23      topK  topP
    14202.0000  5393.6062  0.0000  0.5087  0.0860  0.0000  0.0269  0.0000  0.0136  0.0000  0.0000  0.0179  0.0153  0.0000  0.0490  0.0000  0.0982  0.0000  0.0911  0.0199  0.0143  0.0189  0.0373  0.0000  0.0030  0.0000  1     0.5087
    ```
    * `x` and `y` (float): X Y coordinates of the hexagon.
    * Columns `0` to `{n_factor-1}` (float): Posterior probabilities for each factor `(P(factor | hexagon))`
    * `topK` (int): Index of the most probable factor
    * `topP` (float): Posterior probability of the top factor
* `t{width}_f{n_factor}.bulk_chisq.tsv`
    ```text
    gene  factor  Chi2      pval      FoldChange  gene_total  log10pval
    App   0       864182.5  0.00e+00  4.41        1422353     187657.91
    ```
    * `gene` (str): Gene names.
    * `factor` (int): Factor IDs.
    * `Chi2` (float): Chi-squared test statistic comparing the expression in the target factor and the rest.
    * `pval` (float): P-value associated with the chi-squared test.
    * `FoldChange` (float): Ratio of gene expression inside versus outside the factor.
    * `gene_total` (int): Total count of the gene in the dataset.
    * `log10pval` (float): -log10(pval), useful for ranking significant genes.

* `t{width}_f{n_factor}.factor.info.tsv` and `t{width}_f{n_factor}.factor.info.html`: A TSV and HTML file providing the information for each factor per line.
    ```text
    Factor  RGB          Weight   PostUMI   TopGene_pval                                                                                                                                       TopGene_fc                                                                                                                                           TopGene_weight
    0       255,101,101  0.11279  82949360  App, Snrpn, Calm3, Eef1a2, Dctn1, Ptprn2, Akt3, Tmem181a, Sgtb, Ncdn, Selenow, Atp6v0a1, Cacng8, Atp2a2, Stub1, Akap5, Cmip, Cap1, Rundc3a, Ptprn  Sgtb, Map3k9, Cmip, Vps51, Ptprn2, Zscan2, Ppp2cb, Cacng8, Als2, Eif4e3, Ubxn2a, Gabrb2, Pms1, Rnpc3, Prepl, Zfp418, Tmem181a, Akap5, Dctn1, Atp2b3  App, Calm3, Selenow, Camk2n1, Snrpn, Ncdn, Atp2a2, Eef1a2, Uchl1, Ndfip1, Snap25, Stxbp1, Ywhag, Arf3, Serinc1, Atp6v1b2, Maged1, Atp6v0a1, Dnm1, Rab6b
    ```
    * `Factor` (int): Factor IDs
    * `RGB` (str): Comma-separated RGB values
    * `Weight` (float): Proportion of the total factor signal explained by this factor
    * `PostUMI` (int): Sum of posterior UMI counts across all spatial locations for this factor
    * `TopGene_pval`, `TopGene_fc`, `TopGene_weight` (str): Top marker genes per factor ranked by significance (p-value), fold change, or weight
* `t{width}_f{n_factor}.cmap.tsv`: The color map (the first `n_factor` rows of `--cmap-file`).
    ```text
    R    G    B    Color_hex  Name
    255  101  101  #ff6565    0
    ```
    * `Name` (int): Factor IDs
    * `Color_hex` (str): Color HEX code
    * `R`, `G`, `B` (int): Red, Green, and Blue channel values (0–255)
* `t{width}_f{n_factor}.done` and `t{width}_f{n_factor}_summary.done`: done flags of training and of the DE/report step.

### Decoding Output

The decode ID is `{model_id}_pixel` (e.g., `t18_f24_pixel`).

* `{model_id}_pixel.bin` and `{model_id}_pixel.index`: Binary pixel-level decoding results (top factors and their posterior probabilities per pixel) and their index, read by `punkst draw-pixel-factors` and by packaging.
* `{model_id}_pixel.pseudobulk.tsv.gz`: A feature-by-factor matrix of posterior counts.
    ```text
    Feature  0         1       2       3       4       5       6          7       8       9       10         11        12      13      14        15      16      17      18       19      20       21      22      23
    Neu4     648.5156  0.0000  0.9999  0.0000  0.0000  2.1425  1316.2651  0.0000  0.0000  0.0000  1114.2935  527.5519  0.0000  0.0000  820.3735  0.0000  0.0000  3.3504  30.4680  0.0000  38.0400  0.0000  0.0000  0.0000
    ```
    * Rows (str): Features (e.g., gene names)
    * Columns (int): Factor IDs
    * Values (float): Posterior count of a feature in a factor
* `{model_id}_pixel.bulk_chisq.tsv`: Same format as `t{width}_f{n_factor}.bulk_chisq.tsv`.
* `{model_id}_pixel.factor.info.tsv` and `{model_id}_pixel.factor.info.html`: Same format as `t{width}_f{n_factor}.factor.info.tsv`/`.html`.
* `{model_id}_pixel.png`: An image of the pixel-level factors.
* `{model_id}_pixel.done` and `{model_id}_pixel_summary.done`: done flags.

!!! warning "Use an `--out-dir` path without dots"
    The decode output prefix is derived by `punkst` from the last `.` in the full output path, so an `--out-dir` containing a `.` (e.g., `./out`, `../out`, `run.v2/`) misplaces the decode outputs and the step fails. Use a path without dots (an absolute path such as `/data/run_v2/ficture`).

### UMAP Output

* `t{width}_f{n_factor}.umap.tsv.gz`: UMAP coordinates of the hexagons with their factor assignments; reused by packaging.
* `t{width}_f{n_factor}.umap.png`: A UMAP plot colored by factor.
* `t{width}_f{n_factor}.umap.single.prob.png`: A gallery of UMAP plots, each showing one factor colored by probability.

### 10x Segmentation Output

* `hexagon.d_{width}.10x/`: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`, and the done flag `hexagon_10x.done`.

### Summarization Output

* `ficture.params.json`: A JSON file summarizing input paths and models, read by [`run_cartload2`](./run_cartload2.md).
    ```json
    {
        "in_sge": {
            "in_tiled": "/path/to/out/transcripts.tiled",
            "in_feature": "/path/to/out/transcripts.tiled.selected_features.tsv",
            "in_minmax": "/path/to/sge/coordinate_minmax.tsv"
        },
        "train_params": [
            {
                "model_type": "lda",
                "model_id": "t18_f24",
                "model_path": "/path/to/out/t18_f24.model.tsv",
                "train_width": 18,
                "n_factor": 24,
                "cmap": "/path/to/out/t18_f24.cmap.tsv",
                "decode_params": [
                    {
                        "decode_id": "t18_f24_pixel",
                        "fit_width": 18,
                        "anchor_res": 6,
                        "pixel_res": 0.5,
                        "decode_scale": 1.0
                    }
                ],
                "feature": "/path/to/out/transcripts.tiled.selected_features.hdr.tsv",
                "umap": {
                    "tsv": "/path/to/out/t18_f24.umap.tsv.gz",
                    "png": "/path/to/out/t18_f24.umap.png",
                    "ind_png": "/path/to/out/t18_f24.umap.single.prob.png"
                }
            }
        ]
    }
    ```
    * `in_sge`: the tiled transcript prefix, the FICTURE feature list, and the coordinate min/max file.
    * `train_params`: one entry per model; `feature` is the feature file the model used, `decode_params` lists the decodes (`pixel_res` is `0` with `--single-molecule`), and `umap` appears with `--umap`.

!!! warning "Pass `--in-minmax` (or `--in-json`) with `--tile`"
    When `--tile` runs and `--in-minmax` is omitted, `in_minmax` in `ficture.params.json` is currently recorded as the `transcripts.tiled.done` flag file instead of `transcripts.tiled.coord_range.tsv`. Giving `--in-minmax` (or `--in-json`) avoids this.
