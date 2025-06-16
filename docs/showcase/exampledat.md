# A Quick Start

## Overview

This page provides a quick-start example with a small dataset to help users verify their installation and get familiar with `cartloader` and the `cartostore` workflow.

The workflow consists of four major stages:

1.	SGE Harmonization — prepares raw SGE into a harmonized format.
2.	FICTURE Analysis — computes factor embeddings and visualizations via `punkst`.
3.	Cartloader Compilation — generates web-compatible tiles using `tippecanoe` and `pmtiles`.
4.	AWS Upload — places results into an AWS S3 bucket for sharing via `cartostore`.

## Input Data

!!! info "Input Data Requirement:"
    Input data should be a transcript-indexed SGE containing at least:

    * Spatial coordinates (X coordinates, Y coordinates)
    * Feature metadata (gene symbols)
    * Expression Counts

In this guide, we use an example SGE derived from mouse hippocampus (`Input.tsv.gz`):

*Data Source*: The original data, generated from a full coronal brain slice via SeqScope, was spatially masked to extract a subregion contains hippocampus, returning a compressed TSV (`input.tsv.gz`).

*Data format*:
  ```text 
  #lane  tile  X        Y        gene_id             gene     gn  gt  spl  unspl  ambig
  1      1     7800145  3509475  ENSMUSG00000067786  Nnat     1   1   1    0      0
  1      1     7800145  3510101  ENSMUSG00000020315  Sptbn1   1   1   1    0      0
  1      1     7800145  3520153  ENSMUSG00000024845  Tmem134  1   1   1    0      0
  ```

  * `lane` and `tile`: the lane and tile information
  * `X` and `Y`: Spatial coordinates (pixels, depending on units-per-um)
  * `gene_id`: Ensembl or platform-specific gene identifier (optional)
  * `gene`: Gene Name
  * `gn`, `gt`, `spl`, `unspl`, and `ambig`: Count per gene per pixel for each genomic feature, in the order of Gene, GeneFull, Spliced, Unspliced, and Ambiguous.
      * `gn`: represents unique, confidently mapped transcript count ("gene name"-based);
      * `gt`: denotes total transcript count assigned to gene (includes ambiguities).

*Data Access*:
The input example data can be dowloaded via AWS? Zenodo?

<!-- TODO: decide where to host the input file which is more than 50Mb -->
```bash

```

## Process
### Set Up the Environment

```text
# Define paths to required binaries and resources
spatula=/path/to/spatula/binary         # path to spatula executable
punkst=/path/to/punkst/binary           # path to FICTURE2/punkst executable
tippecanoe=/path/to/tippecanoe/binary   # path to tippecanoe executable
pmtiles=/path/to/pmtiles/binary         # path to pmtiles executable
aws=/path/to/aws/cli/binary             # path to AWS CLI binary

# (Optional) Define path to color map. 
cmap=/path/to/color/map                 # Path to the fixed color map used for rendering. cartloader provides a fixed color map at cartloader/assets/fixed_color_map_256.tsv.

# AWS S3 target location for cartostore
AWS_BUCKET="EXAMPLE_AWS_BUCKET"         # replace this with your actual S3 bucket name

# Unique identifier for your dataset
DATA_ID="seqscope_hippo"                # change this to reflect your dataset name
```

### SGE Harmonization

Convert the raw input to the unified SGE format. See more details in [SGE Harmonization](../step_by_step/sge_harmonization.md).


```bash
cartloader sge_convert \
    --makefn sge_convert.mk \           # (optional) file name of the output make file
    --platform generic \                # use the 'generic' platform parser (adapt as needed for others like 10x_visium_hd, seqscope etc.)
    --in-csv ./input.tsv.gz \           # path to the input.tsv.gz containing raw transcript-indexed SGE
    --csv-colnames-count gn \           # column name for expression counts in the input file (use 'gn' for unique counts in the example data)
    --csv-colname-feature-name gene \   # column name for gene symbols in the input file
    --units-per-um 1000.0 \             # scale to convert coordinates to microns (the example input data is in nanometers, use 1000.0 since 1000 nm = 1 µm)
    --out-dir ./sge \                   # path to output directory where the unified SGE will be saved
    --colnames-count count  \           # output column name for expression count
    --sge-visual \                      # (optional) enable SGE visualization step
    --spatula ${spatula} \              # (optional) path to the spatula binary
    --n-jobs 10                         # (optional) number of parallel jobs for processing
```

### `FICTURE` analysis

Compute spatial factors using `punkst` (FICTURE2 mode). See more details in [FICTURE Analysis](../step_by_step/run_ficture2.md).
```bash
cartloader run_ficture2 \
    --makefn run_ficture2.mk \                          # (optional) file name of the output make file
    --main \                                            # run all five steps in `run_ficture2`
    --in-transcript ./sge/transcripts.unsorted.tsv.gz \ # path to input transcript-level SGE file
    --in-feature ./sge/feature.clean.tsv.gz \           # (optional) path to input feature file
    --in-minmax ./sge/coordinate_minmax.tsv \           # (optional) path to input minmax file
    --cmap-file ${cmap} \                               # (optional) path to input color map file
    --colname-count count \                             # column name for expression count in the input transcript-level SGE file
    --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \        # regex pattern to exclude features (removing mitochondrial and predicted genes in the example analysis)
    --out-dir ./ficture2 \                              # path to output directory
    --width 18 \                                        # LDA training hexagon width (comma-separated if multiple widths are applied)
    --n-factor 6,12,18 \                                # number of factors in LDA training (comma-separated if multiple n-factor are applied)
    --spatula ${spatula} \                              # (optional) path to the spatula binary
    --ficture2 ${punkst} \                              # (optional) path to the punkst directory
    --n-jobs 10  \                                      # (optional) number of parallel jobs 
    --threads 10                                        # (optional) number of threads per job
```

### `cartloader` Compilation
Generate pmtiles and web-compatible tile directories. See more details in [cartloader Compilation](../step_by_step/run_cartload2.md).

```bash
cartloader run_cartload2 \
    --makefn run_cartload2.mk \         # (optional) file name of the output make file
    --fic-dir ./ficture2 \              # path to input directory containing FICTURE2 output
    --out-dir ./cartload2 \             # path to output directory for PMTiles and web tiles
    --id ${DATA_ID} \                   # dataset ID used for naming outputs and metadata
    --colname-count count \             # column name for expression count in the input FICTURE2 output
    --spatula ${spatula} \              # (optional) path to the spatula binary
    --pmtiles ${pmtiles} \              # (optional) path to the pmtiles binary
    --tippecanoe ${tippecanoe} \        # (optional) path to the tippecanoe binary
    --n-jobs 10 \                       # (optional) number of parallel jobs
    --threads 10                        # (optional) number of threads per job
```

### Upload to AWS
Copy the generated cartloader outputs to your designated AWS S3 catalog path:

```bash
cartloader upload_aws_by_catalog \
    --in-dir ./cartload2 \                      # path to the input directory containing the cartloader compilation output
    --s3-dir "s3://${AWS_BUCKET}/${DATA_ID}" \  # path to s3 directory hosting those files
    --aws ${aws} \                              # (optional) path to the aws binary
    --n-jobs 10                                 # (optional) number of parallel jobs
```