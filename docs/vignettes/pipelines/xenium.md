# Xenium Pipeline

## Overview

This tutorial walks through end‑to‑end processing of 10x Xenium data with `CartLoader`: converting inputs, running FICTURE, importing cell results and histology, packaging assets, and uploading to AWS for sharing.

## Input Data

---
**Data Structure and Format**

See details of the Xenium Ranger output at [Xenium Ranger Official Documents](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs)

=== "RAW Data Structure"
    ```
    .
    ├── analysis
    │   ├── clustering
    │   │   ├── gene_expression_graphclust
    │   │   │   └── clusters.csv
    │   │   └── ...
    │   ├── diffexp
    │   │   ├── gene_expression_graphclust
    │   │   │   └── differential_expression.csv
    │   │   └── ...
    │   ├── pca
    │   │   └── gene_expression_10_components
    │   │       ├── projection.csv
    │   │       ├── variance.csv
    │   │       └── ...
    │   └── umap
    │       └── gene_expression_2_components
    │           └── projection.csv
    ├── cell_boundaries.csv.gz
    ├── cell_feature_matrix
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── cells.csv.gz
    ├── morphology_focus
    │   ├── morphology_focus_0000.ome.tif
    │   ├── morphology_focus_0001.ome.tif
    │   ├── morphology_focus_0002.ome.tif
    │   └── morphology_focus_0003.ome.tif
    ├── morphology.ome.tif
    ├── nucleus_boundaries.csv.gz
    ├── transcripts.parquet
    └── ...
    ```

=== "SGE Format"
    !!! info "Naming Convention: `transcripts.csv.gz` or `transcripts.parquet`"
      {%
        include-markdown "../../../includes/includemd_vigenettes_inputformat_xenium.md" start="<!--section1-start-->" end="<!--section1-end-->" 
      %}

=== "Cell Analysis Result Format"
    !!! info "Cells"
        Naming Convention: `cells.csv.gz`
        ```
        "cell_id","x_centroid","y_centroid","transcript_counts","control_probe_counts","genomic_control_counts","control_codeword_counts","unassigned_codeword_counts","deprecated_codeword_counts","total_counts","cell_area","nucleus_area","nucleus_count","segmentation_method"
        "aaaagkdm-1",170.85508728027344,2017.2412109375,1,0,0,0,0,0,1,46.285157930105925,NaN,0,"Segmented by boundary stain (ATP1A1+CD45+E-Cadherin)"
        "aaaamcnn-1",141.60569763183594,2481.442138671875,341,0,0,0,0,1,342,111.5359415486455,50.484689332544804,1,"Segmented by boundary stain (ATP1A1+CD45+E-Cadherin)"
        ```

    !!! info "Cell Boundaries"
        Naming Convention: `cell_boundaries.csv.gz`
        ```
        "cell_id","vertex_x","vertex_y","label_id"
        "aaaagkdm-1",169.3625,2013.0126,1
        "aaaagkdm-1",168.08751,2013.65,1
        ```
    !!! info "Clusters"
        Naming Convention: `clusters.csv`
        ```
        Barcode,Cluster
        aaaagkdm-1,1
        aaaamcnn-1,10
        ```
    !!! info "Differentially Expressed Genes"
        Naming Convention: `differential_expression.csv`
        ```
        Feature ID,Feature Name,Cluster 1 Mean Counts,Cluster 1 Log2 fold change,Cluster 1 Adjusted p value,Cluster 2 Mean Counts,Cluster 2 Log2 fold change,Cluster 2 Adjusted p value,Cluster 3 Mean Counts,Cluster 3 Log2 fold change,Cluster 3 Adjusted p value,Cluster 4 Mean Counts,Cluster 4 Log2 fold change,Cluster 4 Adjusted p value,Cluster 5 Mean Counts,Cluster 5 Log2 fold change,Cluster 5 Adjusted p value,Cluster 6 Mean Counts,Cluster 6 Log2 fold change,Cluster 6 Adjusted p value,Cluster 7 Mean Counts,Cluster 7 Log2 fold change,Cluster 7 Adjusted p value,Cluster 8 Mean Counts,Cluster 8 Log2 fold change,Cluster 8 Adjusted p value,Cluster 9 Mean Counts,Cluster 9 Log2 fold change,Cluster 9 Adjusted p value,Cluster 10 Mean Counts,Cluster 10 Log2 fold change,Cluster 10 Adjusted p value,Cluster 11 Mean Counts,Cluster 11 Log2 fold change,Cluster 11 Adjusted p value,Cluster 12 Mean Counts,Cluster 12 Log2 fold change,Cluster 12 Adjusted p value,Cluster 13 Mean Counts,Cluster 13 Log2 fold change,Cluster 13 Adjusted p value,Cluster 14 Mean Counts,Cluster 14 Log2 fold change,Cluster 14 Adjusted p value,Cluster 15 Mean Counts,Cluster 15 Log2 fold change,Cluster 15 Adjusted p value,Cluster 16 Mean Counts,Cluster 16 Log2 fold change,Cluster 16 Adjusted p value,Cluster 17 Mean Counts,Cluster 17 Log2 fold change,Cluster 17 Adjusted p value,Cluster 18 Mean Counts,Cluster 18 Log2 fold change,Cluster 18 Adjusted p value,Cluster 19 Mean Counts,Cluster 19 Log2 fold change,Cluster 19 Adjusted p value,Cluster 20 Mean Counts,Cluster 20 Log2 fold change,Cluster 20 Adjusted p value,Cluster 21 Mean Counts,Cluster 21 Log2 fold change,Cluster 21 Adjusted p value
        ENSG00000166535,A2ML1,0.008082077992052654,2.292429464470433,0.000013831178435443352,0.0015317620832224837,-0.2717658073791789,0.16427980780868393,0.001226843332106492,-0.42095084864144816,0.46975425102829566,0.0037832334823669654,1.1762824711906763,0.00626099448573822,0.0009370316994608665,-0.9931426193918895,0.00012420457106652693,0.00035171591775058056,-1.9689174885349274,0.002636099946030819,0.0018219664759913663,0.0766263739653219,0.8995675333094153,0.0009475336160580633,-0.8343888826514707,0.06727367639920934,0.004361559793501598,1.5130532082572863,0.0000000000007537933872968394,0.000846294895303063,-0.992227102994244,0.028453359285143767,0.0018146574418219585,0.051194655215844875,0.9304123002179088,0.00045951336258476185,-1.7421234527453535,0.0013628094548944181,0.0023563153416359694,0.4774816523131893,0.2965413340700354,0.0021212060341883753,0.30481838952695206,0.4926601159653632,0.0014237107155417482,-0.14271909128282623,0.8500806231615298,0.0012562334587749313,-0.3507612852540696,0.6384951496662045,0.002355140186915884,0.5936824059089432,0.6436592400308225,0,-1.0564604491458525,0.697077988750366,0.0037541549892459126,1.6545227197901715,0.5114958184195594,0.0016831868337385438,0.9092558204592383,0.9557130750504969,0.004746991923520686,1.8259931272624748,0.4816191969175982
        ENSG00000127837,AAMP,0.02869137687178692,-0.8208377939888187,0.06947481704069494,0.05848830321855647,0.25792687547836657,0.000740291628904809,0.046865415286467996,-0.1268733751597546,0.6859200637377005,0.046074379195969115,-0.1512653200377807,0.6321779498139667,0.06118265802362128,0.30246337282765,0.0011664010874126614,0.04560583066832528,-0.16756035764778598,0.5238325039655534,0.04598894691123035,-0.1611154993878925,0.345954046837975,0.05053512618976337,-0.016703446006141043,0.9871126969288548,0.059006871437084114,0.22791702658442592,0.05960885422623999,0.039698924179670955,-0.37956919506107756,0.02819556072921432,0.04670928255249721,-0.14454178156232,0.25704529967821615,0.04383757479058628,-0.2287656138945824,0.27728049168587215,0.0438060443058687,-0.22809983947157075,0.3116469017790129,0.04846634393266772,-0.08105776515758922,0.7319004973492986,0.029694537781299317,-0.7865793618218326,0.010641440414438465,0.03062069055763895,-0.7473835397757735,0.004183261794659202,0.04811214953271021,-0.08023857624188224,0.9971333031823112,0.036316690185246366,-0.46370796224802824,0.8420547873597443,0.018770774946229564,-1.3137594249747053,0.5160056939276148,0.038713297175986504,-0.34422202322342965,1,0.045887588594033295,-0.11087224163191323,1
        ```

=== "Morphology images"

    Xenium output comes with tissue morphology images in OME‑TIFF, either nuclei‑only (DAPI) or multimodal (DAPI with cell boundary and interior stains). Each file uses a tiled, pyramidal layout (JPEG‑2000, 16‑bit grayscale) with levels from full resolution down to 256×256 for efficient interactive viewing.

    * `morphology.ome.tif`: A 3D Z-stack of the DAPI image.
    * `morphology_focus_0000.ome.tif`: DAPI image
    * `morphology_focus_0001.ome.tif`: Boundary (ATP1A1/E‑Cadherin/CD45)
    * `morphology_focus_0002.ome.tif`: Interior — RNA (18S)
    * `morphology_focus_0003.ome.tif`: Interior — protein (alphaSMA/Vimentin)
    * `morphology_mip.ome.tif`: DAPI maximum intensity projection (MIP) of the Z‑stack.

___

**Data Access**

The example data is downloaded from 10X Data portal.

TO-DO: Code to download the data will be provided later.

___

## Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

Define data ID and analysis parameters:

```bash
# Unique identifier for your dataset
DATA_ID="xenium-v1-humanlung-cancer-ffpe"   # change this to reflect your dataset name
SCALE=1                                     # coordinate to micrometer scaling factor

# LDA parameters
train_width=18                            # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=24                               # define number of factors in LDA training (comma-separated if multiple n-factor are applied)

# s3 path
S3_DIR=/s3/path/to/s3/dir              # Recommend to use DATA_ID as directory name, such as s3://bucket_name/xenium-v1-humanlung-cancer-ffpe
```

!!! info "How to Define Scaling Factors for Xenium"

    The Xenium example data currently used here provides SGE in micrometer units. Use define scaling factor from coordinate to micrometer as 1.

## Run Pipelines

Below is an example of showing running all modules together. You can customize the actions by flags.

In the following example, we only deployed OME_DAPI image. Alternatively, `CartLoader` supports a `--all-images` to deploy all detected image.

```bash
cartloader run_xenium \
  --load-xenium-ranger \
  --sge-convert \
  --run-ficture2 \
  --import-cells \
  --import-images \
  --run-cartload2 \
  --upload-aws \
  --xenium-ranger-dir /path/to/xenium/ranger/output \
  --out-dir /path/to/out/dir \
  --s3-dir ${S3_DIR} \
  --width ${train_width} \
  --n-factor ${n_factor} \
  --id ${DATA_ID} \
  --image-ids OME_DAPI \
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
    `run_xenium` runs multiple `CartLoader` modules; enable and combine actions with flags.


| `CartLoader` Modules                     | Flags in `run_xenium`  | Actions                                                                       | Prerequisites                                                             |
|------------------------------------------|------------------------|-------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| `load_xenium_ranger`                     | `--load-xenium-ranger` | Summarize Xenium Ranger outputs into JSON                                     | Xenium Output files                                                       |
| [`sge_convert`](./sge_convert.md)        | `--sge-convert`        | Convert SGE to CartLoader format; optional density filter and visuals         | Xenium assets JSON (from `load_xenium_ranger`) or transcript CSV/Parquet  |
| [`run_ficture2`](./run_ficture2.md)      | `--run-ficture2`       | FICTURE analysis                                                              | SGE (from `sge_convert`); FICTURE parameters (`--width`, `--n-factor`)    |
| [`import_xenium_cell`](./import_cell.md) | `--import-cells`       | Import cell points, boundaries, cluster, de;                                  | Xenium assets JSON or manual CSVs; also `--cell-id`                       |
| [`import_image`](./import_image.md)      | `--import-images`      | Import background images (OME‑TIFF) → PNG/PMTiles;                            | Xenium assets JSON or `--ome-tifs`; also `--image-ids`  or `--all-images` |
| [`run_cartload2`](./run_cartload2.md)    | `--run-cartload2`      | Package SGE, optional FICTURE/cells/images into PMTiles; write `catalog.yaml` | SGE, optional FICTURE assets or any imported cell/image assets; `--id`    |
| [`upload_aws`](./upload_aws.md)          | `--upload-aws`         | Upload catalog and PMTiles to S3                                              | `catalog.yaml` (from `run_cartload2`); `--s3-bucket`, `--id`              |
| [`upload_zenodo`](./upload_zenodo.md)    | `--upload-zenodo`      | Upload catalog and PMTiles to Zenodo                                          | `catalog.yaml` (from `run_cartload2`); `--zenodo-token`                   |

**Parameter Requirements by Action Flag**

Below are explanations of the parameters used in the example. For the full list, see the [`run_xenium` reference page](../../reference/run_xenium.md).

| Parameter                          | Required when flags                        | Description                                                                           |
|------------------------------------|--------------------------------------------|---------------------------------------------------------------------------------------|
| `--xenium-ranger-dir`              | `--load-xenium-ranger`                     | Xenium Ranger output directory to scan.                                                |
| `--xenium-ranger-assets`           | optional `--load-xenium-ranger`            | Output JSON manifest path (defaults under `--out-dir`).                               |
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
    | `--btf-tifs`                       | `--import-images` (manual)                 | One or more BTF/TIFF paths (relative to `--space-ranger-dir`).                        |
    | `--title`, `--desc`                | optional with `--run-cartload2`            | Human‑readable catalog title/description.                                             | 
    | `--zenodo-token`                   | `--upload-zenodo`                          | Path to file containing Zenodo access token.                                          |-->


## Outputs

### Xenium Assets JSON File

Example: [`includes/xenium_ranger_assets.human_lung_cancer.json`](../../../includes/xenium_ranger_assets.human_lung_cancer.json)

```json
{% include-markdown "../../../includes/xenium_ranger_assets.human_lung_cancer.json" %}
```

### SGE/FICTURE/Cell/Images assets

TO-DO: path to those assets will be provided to serve the output