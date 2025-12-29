# CartLoader


## A Scalable, Harmonized, and Cloud-friendly Ecosystem

Spatial transcriptomics (ST) technologies have revolutionized our ability to map molecular features with remarkable resolution—down to cellular and subcellular levels. The rapid development of sequencing-based platforms (e.g., [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0), [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq), [Pixel-seq](https://www.cell.com/cell/fulltext/S0092-8674(22)01367-8), [10x Visium HD](https://www.10xgenomics.com/platforms/visium)) and imaging-based platforms (e.g., [10x Xenium](https://www.10xgenomics.com/platforms/xenium), [Vizgen MERSCOPE](https://vizgen.com/merscope-ultra/), [CosMx SMI](https://nanostring.com/products/cosmx-spatial-molecular-imager)) has led to an explosion in the production of ST datasets across diverse tissues and species.

However, the rapid diversification of spatial transcriptomics platforms has led to heterogeneous data formats, inconsistent metadata schemas, and platform-specific analysis pipelines. Moreover, existing tools often lack the scalability to visualize gigapixel-scale data without downsampling.

To address this challenge, we introduce a scalable, harmonized, and cloud-friendly CartoScope ecosystem for spatial omics data across platforms.

## Core Modules 

CartoScope ecosystem is composed of three modules:

!!! info "[**CartLoader**](https://seqscope.github.io/cartloader) – a harmonized, scalable pipeline for data processing."

    This document introduces `CartLoader`, the core engine of our ecosystem. It provides a modular, reproducible tool to harmonize, integrate, analyze, and visualize raw high-resolution ST data across platforms. Usage examples are provided in the [Vignettes](./vignettes/intro.md). A detailed description (action, inputs, outputs, and parameters) is provided for each module in the [Reference Pages](./reference/intro.md).

    **Key Features**

    - **Cross-Platform Format Conversion**: Converts raw SGE data from diverse ST platforms into a unified format, enabling consistent downstream processing.
    - **Spatial Factor Inference**: Applies [FICTURE](https://www.nature.com/articles/s41592-024-02415-2) to infer spatial factors directly from pixel-level data, capturing biological patterns without requiring cell segmentation.
    - **Multi-Modal Alignment**: Aligns and overlays histology images with SGE data so all layers share a common coordinate system for pixel-accurate comparisons.
    - **Cloud-Friendly Outputs**: Produces compact, geospatially-indexed data formats suitable for web visualization and cloud storage.
    - **Batch Integration and Sample Stitching**: Supports joint analysis and SGE stitching across samples or platforms to reveal shared or differential features across tissues.
    - **Modular and Reproducible Workflow**: Orchestrates all steps through a Makefile-based system to ensure scalability, transparency, and reproducibility.

??? info "[**CartoStore**](https://github.com/seqscope/cartostore) – an open-access, cloud-hosted repository for dataset sharing."
   
    **CartoStore** acts as the data repository, hosting datasets processed by CartLoader and making them accessible to CartoScope. 
    
    Core Features:

    - **Ready-to-use Datasets**: Access processed, high-resolution datasets of which have FICTURE analysis results, morphology images, cell analysis results, etc., immediately for interactive exploration.
    - **Broad Coverage**: Spans diverse species, tissues, and disease models for a wide range of research needs.
    - **Growing Repository**: Continuously updated with new datasets as they are processed and added.

    Visit the [CartoStore Repository](https://github.com/seqscope/cartostore) to access public datasets.

??? info "[**CartoScope**](index.md) – a web-based App for interactive exploration at original resolution."

    **CartoScope** is the frontend of our ecosystem, enabling interactive, pixel-level exploration of spatial omics data without requiring programming expertise. 

    Core Capabilities:

    - **Ultra-high-resolution data at scale**: Seamless web exploration for submicron-resolution spatial transcriptomics data with hundreds of millions of spatial pixels, tens of thousands of genes, and billions of transcripts.
    - **Molecular-level Inference with FICTURE**: Powered by [FICTURE](https://www.nature.com/articles/s41592-024-02415-2) analysis to preserve spatial accuracy and complexity, breaking the barrier of histology-based segmentation bias.
    - **Dynamic Multi-layer Exploration**: Interactive multi-layer maps with aligned multimodal datasets—including morphology images, cell-level analysis, molecular inference, and spot summaries—enabling exploration with rich, end-to-end spatial context.
    - **Seamless Integration with Morphology Images**: Support for various morphology images (H&E, DAPI, fluorescence reporters) in diverse formats, seamlessly integrated with omics data.
    - **Interactive ROI Analysis**: Interactively define regions of interest, extract spatial features, and run differential expression analysis.
    - **Save & Share Workspace**: One-click save for your workspace (layers, styling, ROIs) and instant sharing with collaborators via unique links.
    - **Open Access**: Upload, screen, and share datasets worldwide at no cost.

    Visit the [CartoScope](placeholder) to learn more details.

Together, this system provides a unified solution for working with raw pixel-level spatial omics data—without sacrificing resolution or reproducibility.

---

## Citations

* FICTURE: [doi.org/10.1038/s41592-024-02415-2](https://www.nature.com/articles/s41592-024-02415-2)
