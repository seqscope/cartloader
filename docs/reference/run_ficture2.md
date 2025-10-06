# Spatial Factor Inference Analysis using FICTURE

## Overview
Following format conversion, `CartLoader` provides the **`run_ficture2`** module to run spatial factor inference using [**`FICTURE`**](https://www.nature.com/articles/s41592-024-02415-2) (Si et al., *Nature Methods*, 2024). This method infers spatial factors directly at the pixel level with submicron resolution, eliminating the need for segmentation.

!!! info "**What is `FICTURE`?**"

    [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2) reconstructs the fine-scale tissue structure by first decomposing gene expression patterns across the tissue section into spatial factors and then assigns each pixel to these factors using local context. Biologically, these inferred factors may correspond to specific cell types, functional or physiological states, subcellular domains, or extracellular transcriptomic signatures.

    By default, [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2) learns the spatial factors by implementing a standard latent Dirichlet allocation (LDA) model on a hexagonal grid overlay of the spatial coordinates. Optionally, spatial factors can also be derived from external sources, such as single-cell or single-nucleus RNA-seq reference datasets, or from spatially agnostic factor learning methods (e.g., Seurat, Scanpy).


!!! info "**The `punkst` version of `FICTURE`**"

    To efficiently run FICTURE-based inference, `CartLoader` integrates [`punkst`](https://github.com/Yichen-Si/punkst), an optimized implementation of FICTURE that maintains output equivalence while enhancing computational scalability and performance. 
    
    Currently, `run_ficture2` uses the `punkst` version of `FICTURE`.

---
## Requirements

- An SGE in the unified format (from [SGE format conversion](./sge_convert.md))
- Pre-installed tools: `spatula`, `punkst`, `gzip`, `sort`, `python`

---
## Example Usage

```bash
cartloader run_ficture2 \
    --main \
    --in-transcript /path/to/converted/transcripts/tsv/file \
    --in-feature /path/to/converted/feature/tsv/file \
    --in-minmax /path/to/converted/coordinates/minmax/tsv/file \
    --cmap-file /path/to/cartloader/assets/fixed_color_map_256.tsv \
    --colname-count count \
    --out-dir /path/to/output/directory \
    --width 18 \
    --n-factor 24 \
    --spatula /path/to/spatula/binary \
    --ficture2 /path/to/punkst/directory \  
    --exclude-feature-regex '^(mt-|Gm\d+$)' \
    --n-jobs 20  \
    --threads 20 
```

---
## Actions

!!! warning "Action Specifications"
    No action runs by default. Activate at least one using the [action parameters](#action-parameters).

### Tiling Step (`--tile`)
The tiling step takes the standardized SGE (from [SGE format conversion](./sge_convert.md)) as input. It reorganizes input coordinate data into non‑overlapping square tiles in a plain TSV format and generates an index file with tile offsets to enable efficient random access.

### Segmentation Step (`--segment`)
The segmentation step starts from the tiled SGE, using the plain TSV file from [tiling step](#tiling-step---tile) as input. It aggregates tiled pixel data into non-overlapping hexagons in a TSV file for spot-level analysis, outputting a tab-delimited file of hexagon records and associated metadata in JSON format.

### LDA Training Step (`--init-lda`)
The LDA training step uses the hexagon TSV and JSON file from [segmentation step](#segmentation-step---segment) as input, trains a Latent Dirichlet Allocation (LDA) model on sparse gene count data from hexagon units, using metadata to interpret input structure and optionally filter or weight features, producing a factorized topic model in TSV.

### Decoding Step (`--decode`)
The decoding step applies a trained LDA model from [LDA training step](#lda-training-step---init-lda) to tiled pixel-level transcript data from [tiling step](#tiling-step---tile) to infer the top spatial factors and their posterior probabilities for each pixel, enabling fine-grained spatial mapping of gene expression. It outputs a pixel-level annotation file in TSV format with coordinates and factor assignments, along with a pseudobulk gene-by-factor matrix in TSV format.

### UMAP Visualization Step (`--umap`)

The UMAP visualization step generate UMAP plots to visualize the relationships between different spatial factors from [LDA training step](#lda-training-step---init-lda).

---
## Parameters

Below are the core parameters. See more details in the collapsible section ("Auxiliary Paramaters") below.

#### Action Parameters


* `--main`: Run all of the following five actions.
* `--tile`: Run [tiling step](#tiling-step---tile).
* `--segment`: Run [segmentation step](#segmentation-step---segment).
* `--init-lda`: Run [LDA training step](#lda-training-step---init-lda).
* `--decode`: Run [decoding step](#decoding-step---decode).
* `--umap`: Run [UMAP step](#umap-visualization-step---umap).

#### Input/Output Parameters

* `--out-dir` (str): Output directory to store all result files.
* `--out-json` (str): Output JSON file summarizing FICTURE parameters (Default: `<out-dir>/ficture.params.json`).
* `--in-transcript` (str): Input transcript-indexed SGE file in TSV format.
* `--in-minmax` (str): Optional input coordinate min-max file.
* `--in-feature` (str): Optional input UMI count per gene TSV file.

#### Key Parameters

* `--width` (str): Comma-separated hexagon flat-to-flat widths (in µm) for LDA training.
* `--n-factor` (str): Comma-separated list of factor counts for LDA training.
* `--include-feature-regex` (str): (Optional) Regex pattern for including features/genes.
* `--exclude-feature-regex` (str): (Optional) Regex pattern for excluding features/genes.
* `--cmap-file` (str): (Optional) Path to fixed color map TSV file. If not provided, FICTURE will generate a color map.

??? note "Auxiliary Paramaters"

    **Auxiliary Input Parameters**

    * `--in-feature-ficture` (str): (Optional) A separate feature file applied to FICTURE analysis. Alternative to customizing via auxiliary parameters.
    * `--colidx-x` (int): Column index of X in the transcript file (Default: 1).
    * `--colidx-y` (int): Column index of Y in the transcript file (Default: 2).
    * `--colname-count` (str): Column name to use as count value (Default: count).
    * `--colname-feature` (str): Column name for gene/feature name (Default: gene).

    **Auxiliary FICTURE Parameters**:

    * Tiling-specific parameters:
        * `--tile-size` (int): Size of tiles for processing (Default: 500).
        * `--tile-buffer` (int): Buffer zone around each tile (Default: 1000).
    * Segmentation-specific parameters:
        * `--min-ct-per-unit-hexagon` (int): Minimum count per hexagon (Default: 50).
        * `--minibatch-size` (int): Minibatch size for preprocessing (Default: 500).
    * LDA training-specific parameters:
        * `--min-ct-per-unit-train` (int): Minimum count for training (Default: 50).
        * `--train-epoch` (int): Number of epochs to train LDA model (Default: 2).
    * Decoding-specific parameters:
        * `--fit-width` (int): Hexagon width (in µm) for model fitting (Default: same as train width).
        * `--min-ct-per-unit-fit` (int): Minimum count per unit during model fitting (Default: 50).
        * `--anchor-res` (int): Anchor resolution used in decoding (Default: 6).
        * `--radius-buffer` (int): Buffer added to anchor resolution for decoding (Default: 1).
        * `--fit-plot-um-per-pixel` (int): Image resolution for fit coarse plots (Default: 1).
        * `--decode-scale` (int): Scales input coordinates to pixels in the output image (Default: 1)')
    * Shared paramaters across steps:
        * `--seed` (int): Random seed for reproducibility (Default: 1).
        * `--min-ct-per-feature` (int): Minimum count per feature for LDA and decoding (Default: 20).
        * `--de-max-pval` (float): p-value cutoff for differential expression (Default: 1e-3).
        * `--de-min-fold` (float): Fold-change threshold for differential expression (Default: 1.5).

    ** Environment Parameters**:
    For tools that require specifying the path to their executable binaries, you may omit the path if the binary is already included in your system's `PATH`.

    * `--gzip` (str): Path to `gzip` binary; consider `pigz -p 4` for speed (Default: `gzip`).
    * `--sort` (str): Path to `sort` binary (Default: `sort`).
    * `--sort-mem` (str): Memory allocated per `sort` process (Default: 1G).
    * `--spatula` (str): Path to `spatula` binary (Default: `spatula`).
    * `--ficture2` (str): Path to the `punkst` repository (Default to `punkst` directory in `submodules`)
    * `--python` (str): Path to Python 3 binary (Default: `python3`).

    ** Run Parameters**:

    * `--dry-run` (flag): Generate the Makefile; do not execute.
    * `--restart` (flag): Ignore existing outputs and rerun all steps.
    * `--makefn` (str): Makefile name to write (default: `run_ficture2.mk`).
    * `--n-jobs` (int): Number of parallel jobs (default: 1).
    * `--threads` (int): Max threads per job (default: 4).

---
## Output

### Tiling Output

* `transcripts.tiled.tsv`: A tab-delimited TSV file where each line records a tile with X and Y coordinates, a feature name (e.g., gene), and its associated count.
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
* `transcripts.tiled.features.tsv`: A tab-delimited TSV file listing each feature name and its expression count per line.
    ```text
    1810009J06Rik	1
    Gm48545	        1
    ```
### Segmentation Output

* `hexagon_d_{width}.randomized.tsv`: A plain-text TSV file that stores sparse feature representations for hexagon units, with each line corresponding to one hexagon.
    
    !!! info "File Format"
        This is not a regular table—the number of columns varies across lines depending on the number of features present in each hexagon.

    ```text
    00001c64  18999.0000  9773.9627  227  244   2417 1  10469 2  2448 1  543 1   10935 1  6460 1  1797 1  6517 1  1094 1  9872 1  8587 1  7137 1  8786 1  10564 1  2127 1  439 1   1529 1  8432 1   1304 1  1739 1   3012 1  6434 1  5297 1  10017 1  2567 1  4771 1  5644 1  6415 1  5258 1   258 1   2314 1  6378 1  8843 1  7386 1  3322 1  2041 1  4735 1  9755 1  8614 1   4665 1  5693 1  2194 1  1680 1  10406 1  2696 1  7250 1  2367 1  7189 1  1616 1  10721 1  441 1   9981 1   8696 1  1012 1  8545 1  835 1   6550 1  4409 1  1860 1  575 1   139 1   5473 1  10172 1  7316 1  5260 1  180 1   2068 1  10806 1  8817 1  1621 1  4542 1  6193 1  440 1   6608 1  6989 1  4782 1  3258 1  5205 1  1093 1  9060 1  1921 1  6333 1   846 1   10575 1  2608 1  10186 1  932 1   2217 1  5392 1  2275 1  4845 1   8401 1  341 1   2810 1   754 1   4501 1  1160 1  5354 1  5097 1  3713 1  1619 1  9843 1  2749 1  1207 1  7384 1  2758 1  445 1   9040 1  45 1    2517 1   11255 1  8161 1  393 1   2963 1  3457 1  4481 1  10649 1  6954 2  2842 1  7982 1  2034 1  2849 2  5393 1  5257 1  773 1   1030 1  5455 1  2380 1  7079 1  1491 1  353 1    5343 1  1231 1  178 1   8659 1  10994 1  6168 1  267 1    6770 1   9083 1   2691 1  2423 2  3058 3  6921 1  1778 1  236 1   7175 1  7683 1   47 1    9094 2  1898 1  6524 1  5640 1  6668 1   6411 1  3536 1  10475 2  6538 1  1443 1  2728 1  7717 1  4633 1  6634 1   5349 1  11215 2  9159 1  2129 1  5984 1  8584 2  7042 1  3819 2  7931 1  2343 1  1182 2  2724 1  1871 1  716 1   1744 1  263 1   7459 1  6737 1  4654 1   9023 1  775 2    1546 1  1512 1  6395 1  5624 1  1560 2  8919 2  4689 1  2403 1  4163 1  5705 1  1260 1  8068 1  1900 1  1701 1  1187 1  8172 1  5955 1  494 1   428 1   6082 1  479 1   5362 1   11273 1  1764 1   9546 1  8261 1  11088 1  6262 1  2154 1  612 1   6777 2  866 1    10850 1  4168 1  9334 1  9077 2   3818 1  732 1   6386 1  989 1   1277 1  319 1   6230 1
    ```
	
    * 1st column (str): Random identifier for the hexagon
	* 2nd and 3rd columns (float): Axial coordinates in the hexagonal coordinate system
	* 4th and 5th columns (int): Number of unique features and total counts, repeated per modality (for K modalities)
	* 6th column onward (integer integer): Index–count pairs for non-zero features, using 0-based indices from --feature-dict or as provided in the input file

* `hexagon_d_{width}.json`: The meta data for the hexagons. such as, , hexagon size, x and y column indices (), number of features, number of modalities, and number of hexagons (`n_units`),
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
    * `dictionary`: The feature name-idex pairs
    * `header_info`: Header information
    * `hex_size`: Size of hexagons
    * `icol_x` and `icol_y`: Indices of X and Y coordinates
    * `n_features`: Number of features
    * `n_modalities`: Number of modalities
    * `n_units`: Number of hexagons
    * `offset_data`: ?
    * `random_key`: the random seed for reproduction

### LDA Training Output

* `t{width}_f{n_factor}.model.tsv`: A matrix to store the learned topic-word distribution.
    ```text
    Feature  0        1      2      3      4      5      6         7      8      9      10        11       12     13     14        15     16     17     18      19     20      21     22     23
    Neu4     481.168  0.042  0.042  0.042  0.042  0.042  1828.293  0.042  0.042  0.042  2290.443  836.722  0.042  0.042  1704.589  0.042  0.042  0.042  29.738  0.042  24.854  0.042  0.042  0.042
    ```
    * Rows (str): LDA topics or spatial factors
	* Columns (int): Feature identifiers (e.g., gene names or indices)
	* Values (float): Probability or weight of a feature under each topic `(P(feature | topic))`
* `t{width}_f{n_factor}.results.tsv.gz`: A tab-delimited file containing the posterior topic distributions for each hexagon.  
    ```text
    x           y          0       1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21      22      23      topK  topP
    14202.0000  5393.6062  0.0000  0.5087  0.0860  0.0000  0.0269  0.0000  0.0136  0.0000  0.0000  0.0179  0.0153  0.0000  0.0490  0.0000  0.0982  0.0000  0.0911  0.0199  0.0143  0.0189  0.0373  0.0000  0.0030  0.0000  1     0.5087
    ```
    * `x` and `y` (float): X Y coordinates of the hexagon.
    * Columns `0` to `{n_factor-1}` (float): Posterior probabilities for each latent factor `(P(topic | hexagon))`
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
    * `FoldChange` (float): Ratio of gene expression inside versus outside the factor’s high-loading region.
    * `gene_total` (int): Total count of the gene in the dataset.
    * `log10pval` (float): Base-10 logarithm of the inverse p-value (i.e., -log10(pval)), useful for ranking significant genes.

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
* `t{width}_f{n_factor}.cmap.tsv`: a color map.
    ```text
    R    G    B    Color_hex  Name
    255  101  101  #ff6565    0
    ```
    * `Name` (int): Factor IDs
    * `Color_hex` (str): Color HEX code
    * `R`, `G`, `B` (float): Red, Green, and Blue channel values (range: 0.0 to 1.0)

### Decoding Output

* `t{width}_f{n_factor}_p{width}_a{anchor_res}.tsv.gz`: A tab-delimited file where each line represents a pixel–feature pair, recording the pixel’s coordinates, the expressed feature and its count, along with the top three most probable latent factors (`K1`–`K3`) and their corresponding probabilities (`P1`–`P3`).
    ```text
    #x         y         feature  ct  K1  K2  K3  P1          P2          P3
    5159.4800  450.1500  Retreg2  1   0   19  9   9.5073e-01  4.8011e-02  5.0298e-04
    ```
    * `x` and `y` (float): X Y coordinates
    * `feature` (str): feature names.
    * `ct` (int): count
    * `K1` (int) and `P1` (float): The most probable factor and its probability
    * `K2` (int) and `P2` (float): The 2nd most probable factor and its probability
    * `K3` (int) and `P3` (float): The 3rd most probable factor and its probability

* `t{width}_f{n_factor}_p{width}_a{anchor_res}.index`: An index file for `t{width}_f{n_factor}_p{width}_a{anchor_res}.tsv.gz`.
* `t{width}_f{n_factor}_p{width}_a{anchor_res}.json`: A JSON file to provide header information for `t{width}_f{n_factor}_p{width}_a{anchor_res}.tsv.gz`.
    ```json 
    {
        "K1": 4,
        "K2": 5,
        "K3": 6,
        "P1": 7,
        "P2": 8,
        "P3": 9,
        "ct": 3,
        "feature": 2,
        "x": 0,
        "y": 1
    }
    ```
* `t{width}_f{n_factor}_p{width}_a{anchor_res}.png`: A PNG file to visualize the spatial factors distribution at pixel level.
* `t{width}_f{n_factor}_p{width}_a{anchor_res}.pseudobulk.tsv`: A feature-by-factor matrix showing feature distribution across topics.
    ```text
    Feature  0         1       2       3       4       5       6          7       8       9       10         11        12      13      14        15      16      17      18       19      20       21      22      23
    Neu4     648.5156  0.0000  0.9999  0.0000  0.0000  2.1425  1316.2651  0.0000  0.0000  0.0000  1114.2935  527.5519  0.0000  0.0000  820.3735  0.0000  0.0000  3.3504  30.4680  0.0000  38.0400  0.0000  0.0000  0.0000
    ```
    * Rows (str): LDA topics or spatial factors
	* Columns (int): Feature identifiers (e.g., gene names or indices)
	* Values (float): Posterior probability of a feature in a factor.

* `t{width}_f{n_factor}_p{width}_a{anchor_.res}.bulk_chisq.tsv`: Same format as the `t{width}_f{n_factor}.bulk_chisq.tsv`
* `t{width}_f{n_factor}_p{width}_a{anchor_res}.factor.info.tsv`:  Same format as the `t{width}_f{n_factor}.factor.info.tsv`
* `t{width}_f{n_factor}_p{width}_a{anchor_res}.factor.info.html`: Same format as the `t{width}_f{n_factor}.factor.info.html`

### Summarization Output

* `ficture.params.json`: A JSON file to summarize the path to input, output, and paramaters.
    ```json
    {
    "in_sge": {
        "in_transcript": "/path/to/transcripts.unsorted.tsv.gz",
        "in_feature": "/path/to/transcripts.tiled.features.hdr.tsv",
        "in_minmax": "/path/to/coordinate_minmax.tsv"
    },
    "in_feature_ficture": "/path/to/transcripts.tiled.overlapping_features.min100.hdr.tsv", # if available
    "train_params": [
        {
            "model_type": "lda",
            "model_id": "t{width}_f{n_factor}",
            "model_path": "/path/to/t{width}_f{n_factor}.model.tsv",
            "train_width": width,
            "n_factor": n_factor,
            "cmap": "/path/to/t{width}_f{n_factor}.cmap.tsv",
            "decode_params": [
                {
                    "decode_id": "t{width}_f{n_factor}_p{width}_a{anchor_res}",
                    "fit_width": width,
                    "anchor_res": anchor_res
                }
            ]
        }
    ]
    }
    ```

### UMAP Output

* `t{width}_f{n_factor}.umap.png`: A UMAP plot that visualizes the relationships between spatial factors.
* `t{width}_f{n_factor}.umap.single.prob.png`: A gallery of images, with each image showing a single factor color-coded by probability.
