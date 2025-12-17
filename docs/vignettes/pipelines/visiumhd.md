# Visium HD End-to-End Pipeline

## Overview

This tutorial walks through end‑to‑end processing of 10x Visium HD data with `CartLoader`: converting inputs, running FICTURE, importing cell results and histology, packaging assets, and uploading for sharing.

## Input Data
---
**Data Structure and Format**

See Space Ranger output details in the official documentation: [Space Ranger Outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview)

=== "RAW Data Structure"
    ```
    ├── binned_outputs
    │   ├── square_002um
    │   │   ├── filtered_feature_bc_matrix
    │   │   │   ├── barcodes.tsv.gz
    │   │   │   ├── features.tsv.gz
    │   │   │   └── matrix.mtx.gz
    │   │   ├── ...
    │   │   └── spatial
    │   │       ├── aligned_fiducials.jpg
    │   │       ├── aligned_tissue_image.jpg
    │   │       ├── cytassist_image.tiff
    │   │       ├── detected_tissue_image.jpg
    │   │       ├── scalefactors_json.json
    │   │       ├── tissue_hires_image.png
    │   │       ├── tissue_lowres_image.png
    │   │       └── tissue_positions.parquet
    │   ├── square_008um
    │   │   ├── analysis
    │   │   │   ├── clustering
    │   │   │   │   ├── gene_expression_graphclust
    │   │   │   │   │   └── clusters.csv
    │   │   │   │   └── ...
    │   │   │   ├── diffexp
    │   │   │   │   ├── gene_expression_graphclust
    │   │   │   │   │   └── differential_expression.csv
    │   │   │   │   └── ...
    │   │   │   ├── pca
    │   │   │   │   └── gene_expression_10_components
    │   │   │   │       ├── projection.csv
    │   │   │   │       ├── variance.csv
    │   │   │   │       └── ...
    │   │   │   └── umap
    │   │   │       └── gene_expression_2_components
    │   │   │           └── projection.csv
    │   │   ├── filtered_feature_bc_matrix
    │   │   │   └── ...     # Mirrors binned_outputs/square_002um/filtered_feature_bc_matrix
    │   │   ├── ...
    │   │   └── spatial
    │   │       └── ...     # Mirrors binned_outputs/square_002um/spatial
    │   └── square_016um
    │       └── ...         # Mirrors square_008um (analysis, spatial, filtered_feature_bc_matrix)
    ├── segmented_outputs
    │   ├── cell_segmentations.geojson
    │   └── ...             # Mirrors binned_outputs/square_008um (analysis, spatial, filtered_feature_bc_matrix)
    └── ...
    ```

=== "SGE Format"
    {%
      include-markdown "../../../includes/includemd_vigenettes_inputformat_visiumhd.md" 
    %}

=== "Cell Analysis Result Format"
    !!! info "Cell Segmentation Mask"
        Location: `segmented_outputs/cell_segmentations.geojson`
        ```
        "cell_id","x_centroid","y_centroid","transcript_counts","control_probe_counts","genomic_control_counts","control_codeword_counts","unassigned_codeword_counts","deprecated_codeword_counts","total_counts","cell_area","nucleus_area","nucleus_count","segmentation_method"
        "aaaagkdm-1",170.85508728027344,2017.2412109375,1,0,0,0,0,0,1,46.285157930105925,NaN,0,"Segmented by boundary stain (ATP1A1+CD45+E-Cadherin)"
        "aaaamcnn-1",141.60569763183594,2481.442138671875,341,0,0,0,0,1,342,111.5359415486455,50.484689332544804,1,"Segmented by boundary stain (ATP1A1+CD45+E-Cadherin)"
        ```
    !!! info "Cell Feature Matrix"
        Location: `segmented_outputs/filtered_feature_cell_matrix`
        
        This contains `barcodes.tsv.gz`,`features.tsv.gz`, `matrix.mtx.gz`, providing cell id, feature information, and expression count matrix, respectively.

    !!! info "Clusters"
        Location: `segmented_outputs/analysis/clustering/gene_expression_graphclust/clusters.csv`
        ```
        Barcode,Cluster
        cellid_000000001-1,22
        cellid_000000002-1,22
        ```
    !!! info "Differentially Expressed Genes"
        Location: `segmented_outputs/analysis/diffexp/gene_expression_graphclust/differential_expression.csv`
        ```
        Feature ID,Feature Name,Cluster 1 Mean Counts,Cluster 1 Log2 fold change,Cluster 1 Adjusted p value,Cluster 2 Mean Counts,Cluster 2 Log2 fold change,Cluster 2 Adjusted p value,Cluster 3 Mean Counts,Cluster 3 Log2 fold change,Cluster 3 Adjusted p value,Cluster 4 Mean Counts,Cluster 4 Log2 fold change,Cluster 4 Adjusted p value,Cluster 5 Mean Counts,Cluster 5 Log2 fold change,Cluster 5 Adjusted p value,Cluster 6 Mean Counts,Cluster 6 Log2 fold change,Cluster 6 Adjusted p value,Cluster 7 Mean Counts,Cluster 7 Log2 fold change,Cluster 7 Adjusted p value,Cluster 8 Mean Counts,Cluster 8 Log2 fold change,Cluster 8 Adjusted p value,Cluster 9 Mean Counts,Cluster 9 Log2 fold change,Cluster 9 Adjusted p value,Cluster 10 Mean Counts,Cluster 10 Log2 fold change,Cluster 10 Adjusted p value,Cluster 11 Mean Counts,Cluster 11 Log2 fold change,Cluster 11 Adjusted p value,Cluster 12 Mean Counts,Cluster 12 Log2 fold change,Cluster 12 Adjusted p value,Cluster 13 Mean Counts,Cluster 13 Log2 fold change,Cluster 13 Adjusted p value,Cluster 14 Mean Counts,Cluster 14 Log2 fold change,Cluster 14 Adjusted p value,Cluster 15 Mean Counts,Cluster 15 Log2 fold change,Cluster 15 Adjusted p value,Cluster 16 Mean Counts,Cluster 16 Log2 fold change,Cluster 16 Adjusted p value,Cluster 17 Mean Counts,Cluster 17 Log2 fold change,Cluster 17 Adjusted p value,Cluster 18 Mean Counts,Cluster 18 Log2 fold change,Cluster 18 Adjusted p value,Cluster 19 Mean Counts,Cluster 19 Log2 fold change,Cluster 19 Adjusted p value,Cluster 20 Mean Counts,Cluster 20 Log2 fold change,Cluster 20 Adjusted p value,Cluster 21 Mean Counts,Cluster 21 Log2 fold change,Cluster 21 Adjusted p value,Cluster 22 Mean Counts,Cluster 22 Log2 fold change,Cluster 22 Adjusted p value,Cluster 23 Mean Counts,Cluster 23 Log2 fold change,Cluster 23 Adjusted p value,Cluster 24 Mean Counts,Cluster 24 Log2 fold change,Cluster 24 Adjusted p value
        ENSMUSG00000051951,Xkr4,0.0020893375742615837,0.542453177216764,0.16305286911324085,0.0019306072327678393,0.34251290321116556,0.628624629681027,0.0009020306879791494,-0.677062452566048,0.6057595443651351,0.0008575468486789855,-0.53825692400323,1,0.0018835945742411177,0.3347722133550999,1,0.0017232930683496192,0.30174489068131116,1,0.0008765004215231041,-0.6062657172996246,1,0.0006958304617047569,-0.669599363115374,1,0.0021901445815864866,0.5424639235287998,0.6460045298567026,0,-0.513593267969263,1,0.000450562773540615,-0.8809942923847434,1,0.0006856395158652269,-0.2651228613954135,1,0.0008800111512863451,-0.49971125727530463,1,0.0026227703011305675,1.0167399161570962,1,0.002277387540951688,0.7536928538206205,1,0,-1.333742893536682,1,0.0007018464992057343,-0.23097850703201495,1,0.0005602354699616488,-0.560849410646286,1,0.0015894120085537071,0.37398251412889394,1,0,1.526558986825357,1,0,-0.5313075119163866,1,0,4.9659921506021,1,0,1.1981124565354708,1,0,1.1280527529013042,1
        ENSMUSG00000089699,Gm1992,0.00003369899313325135,2.372528175774452,0.49812531917443914,0,1.0054779159335947,1,0,2.907900048155108,1,0,3.8795955908826674,1,0,2.9741824980986316,1,0,3.8867073914024672,1,0,3.4811971239507127,1,0,4.171702890865568,1,0,2.5799386289474633,1,0,5.929350227879466,1,0,4.553633935251982,1,0,5.16950536624131,1,0,3.9181412576105927,1,0,5.104202757407435,1,0,4.56960978938165,1,0,5.109200602312047,1,0,5.203649720604709,1,0,4.873778816990438,1,0,4.791835029014791,1,0,7.969502482674086,1,0,5.911635983932342,1,0,11.408935646450828,1,0,7.641055952384198,1,0,7.570996248750033,1
        ```
___

**Data Access**

Downloaded the ST data from [10x Genomics Dataset portal](https://www.10xgenomics.com/datasets/visium-hd-three-prime-mouse-brain-fresh-frozen).

___

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
# Unique identifier for your dataset
DATA_ID="visiumhd_3prime_mouse_brain"    # change this to reflect your dataset name

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)

# Path to AWS S3 directory
S3_DIR=/s3/path/to/s3/dir                # Recommend to use DATA_ID as directory name, such as s3://bucket-name/visiumhd-3prime-mouse-brain
```

!!! info "How to define scaling for Visium HD?"

    10x Visium HD provides `scalefactors_json.json` (pixel‑to‑µm). `CartLoader` accepts it via `--scale-json` and computes the scaling automatically, so you don’t need to manually specify `--units-per-um`.
    
    Alternatively, provide the scale directly with `--units-per-um`.


## Run Pipelines

The example below runs all modules together. Customize actions with flags.

```bash
cartloader run_visiumhd \
  --load-space-ranger \
  --sge-convert \
  --run-ficture2 \
  --import-cells \
  --import-images \
  --run-cartload2 \
  --upload-aws \
  --space-ranger-dir /path/to/space/ranger/output \
  --out-dir /path/to/out/dir \
  --s3-dir ${S3_DIR} \
  --width ${train_width} \
  --n-factor ${n_factor} \
  --id ${DATA_ID} \
  --spatula ${spatula} \
  --ficture2 ${punkst} \
  --pmtiles ${pmtiles} \
  --tippecanoe ${tippecanoe} \
  --aws ${aws} \
  --n-jobs ${n_jobs} \
  --threads ${n_jobs}
```

**Action Flags to Enable Modules**

!!! warning "Actions"
    `run_visiumhd` runs multiple `CartLoader` modules together; enable and combine actions with flags.

| `CartLoader` Modules                    | Flags in `run_visiumhd` | Actions                                                                       | Prerequisites                                                                 |
|-----------------------------------------|-------------------------|-------------------------------------------------------------------------------|-------------------------------------------------------------------------------|
| `load_space_ranger`                     | `--load-space-ranger`   | Summarize Space Ranger outputs into JSON                                      | Space Ranger Output files                                                     |
| [`sge_convert`](./sge_convert.md)       | `--sge-convert`         | Convert SGE to `CartLoader` format; optional density filter and visuals         | Space Ranger assets JSON (from `load_space_ranger`) or transcript CSV/Parquet |
| [`run_ficture2`](./run_ficture2.md)     | `--run-ficture2`        | FICTURE analysis                                                              | SGE (from `sge_convert`); FICTURE parameters (`--width`, `--n-factor`)        |
| [`import_space_cell`](./import_cell.md) | `--import-cells`        | Import cell points, boundaries, cluster, de;                                  | Space Ranger assets JSON or manual CSVs; also `--cell-id`                     |
| [`import_image`](./import_image.md)     | `--import-images`       | Import background images (BTF-TIFF) → PNG/PMTiles;                            | Space Ranger assets JSON or `--tifs`; also `--image-ids`/`--all-images`   |
| [`run_cartload2`](./run_cartload2.md)   | `--run-cartload2`       | Package SGE, optional FICTURE/cells/images into PMTiles; write `catalog.yaml` | SGE, optional FICTURE assets or any imported cell/image assets; `--id`        |
| [`upload_aws`](./upload_aws.md)         | `--upload-aws`          | Upload catalog and PMTiles to S3                                              | `catalog.yaml` (from `run_cartload2`); `--s3-bucket`, `--id`                  |
| [`upload_zenodo`](./upload_zenodo.md)   | `--upload-zenodo`       | Upload catalog and PMTiles to Zenodo                                          | `catalog.yaml` (from `run_cartload2`); `--zenodo-token`                       |

**Parameter Requirements by Action Flag**

Below are explanations of the parameters used in the example. For the full list, see the [`run_visiumhd` reference page](../../reference/run_visiumhd.md).

| Parameter                          | Required when flags                        | Description                                                                           |
|------------------------------------|--------------------------------------------|---------------------------------------------------------------------------------------|
| `--space-ranger-dir`               | `--load-space-ranger`                      | Space Ranger output directory to scan.                                                |
| `--space-ranger-assets`            | optional `--load-space-ranger`             | Output JSON manifest path (defaults under `--out-dir`).                               |
| `--width`                          | `--run-ficture2`                           | Hexagon width(s) in µm for training/projection.                                       |
| `--n-factor`                       | `--run-ficture2`                           | Factor count(s) for FICTURE training.                                                 |
| `--cell-id`                        | `--import-cells`                           | Asset ID/prefix for cell outputs.                                                     |
| `--image-ids` or `--all-images`    | `--import-images`                          | Choose specific image IDs or import all detected.                                     |
| `--id`                             | `--run-cartload2`                          | Catalog ID (used in catalog.yaml and outputs).                                        |
| `--s3-dir`                         | `--upload-aws`                             | Destination S3 path (e.g., `s3://bucket/prefix`).                                     |
| `--out-dir`                        | any action                                 | Output root directory for generated artifacts.                                        |
| `--dry-run`, `--restart`           | optional (any)                             | Control execution (preview or rerun ignoring existing outputs).                       |
| `--n-jobs`, `--threads`            | optional (any)                             | Parallelism for Make/GDAL/tippecanoe steps.                                           |


<!-- | `--in-json`                        | `--sge-convert` with `--load-space-ranger` | Use `--space-ranger-assets` for SGE inputs (JSON mode).                          |
| `--in-mex`                         | `--sge-convert` (manual)                   | Path to `square_002um/filtered_feature_bc_matrix` (relative to `--space-ranger-dir`). |
| `--pos-parquet`                    | `--sge-convert` (manual)                   | Path to `square_002um/spatial/tissue_positions.parquet` (relative).                   |
| `--scale-json` or `--units-per-um` | `--sge-convert` (manual)                   | Provide pixel-to‑µm scale via `scalefactors_json.json` or numeric factor.             |
| `--filter-by-density`              | optional with `--sge-convert`              | Enable density‑based transcript filtering.                                            |
| `--exclude-feature-regex`          | optional with `--sge-convert`              | Regex to exclude features (e.g., controls).                                           |
| `--csv-clust`, `--csv-diffexp`     | `--import-cells` (manual)                  | Paths to clusters/DE CSVs (relative to `--space-ranger-dir`).                         |
| `--tifs`                       | `--import-images` (manual)                 | One or more BTF/TIFF paths (relative to `--space-ranger-dir`).                        |
| `--title`, `--desc`                | optional with `--run-cartload2`            | Human‑readable catalog title/description.                                             | 
| `--zenodo-token`                   | `--upload-zenodo`                          | Path to file containing Zenodo access token.                                          |-->

## Outputs

<!-- ### Space Ranger Assets JSON

```json
{% include-markdown "../../../includes/visiumhd_space_ranger_assets.visiumhd_3prime_mouse_brain.json" %}
```

### Spatial Factor Inference
Below is an example of spatial factor inference results from `FICTURE` using a training width of 18, 12 factors, a fit width of 18, and an anchor resolution of 6. See more details of output at the Reference pages for [run_ficture2](../docs/reference/run_ficture2.md)

![FICTURE](../../images/pipeline_vignettes/visiumhd_3prime_mouse_brain.t18_f24_p18_a6.png)
![cmap](../../images/pipeline_vignettes/visiumhd_3prime_mouse_brain.t18_f24.rgb.png)

{{ read_csv('../../tabs/visiumhd_3prime_mouse_brain.t18_f24_p18_a6.factor.info.tsv',sep = '\t') }}

### SGE/FICTURE/Cell/Images assets

See more details of output at the Reference pages for [run_cartload2](../..//reference/run_cartload2.md), [import_xenium_cell](../..//reference/import_cell.md), and [import_image](../../reference/import_image.md).

Individual PMTiles and asset JSON files reside alongside it under `<out-dir>/cartload2/`.
- Catalog to serve: `<out-dir>/cartload2/catalog.yaml`
- SGE PMTiles (examples): `<out-dir>/cartload2/sge-mono-dark.pmtiles`, `<out-dir>/cartload2/sge-mono-light.pmtiles`
- Factor PMTiles (examples): `<out-dir>/cartload2/t18_f24_p18_a6-results.pmtiles`, `<out-dir>/cartload2/t18_f24_p18_a6-pixel-raster.pmtiles`
- Cells and Boundaries PMTiles: `<out-dir>/cartload2/spaceranger-cells.pmtiles`, `<out-dir>/cartload2/spaceranger-boundaries.pmtiles` -->

<div class="grid cards generic" markdown>

-   ![Interactive Exploration](../../images/pipeline_vignettes/visiumhd_3prime_mouse_brain.t18_f24_p18_a6.png)

    ---

    #### View/Explore

    The output are available in both CartoScope. 

    [Explore in CartoScope](http://localhost:5173/dataset?uri=s3%2Fcartostore%2Fdata%2Fbatch%3D2025_12%2Fcartloader-pipeline-example-collection%2Fvisiumhd_3prime_mouse_brain){ .md-button .md-button--primary .button-tight-small }

    <!-- [Download from Zenodo](https://zenodo.org/records/17958847){ .md-button .button-tight-small } -->

</div>

See more details of output at the Reference pages for [run_ficture2](../docs/reference/run_ficture2.md) and [run_cartload2](../docs/reference/run_cartload2.md).


