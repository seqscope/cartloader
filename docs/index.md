# A scalable, harmonized, and cloud-friendly ecosystem for spatial transcriptomics data

## Overview

Spatial transcriptomics (ST) technologies have revolutionized our ability to map molecular features with remarkable resolution—down to cellular and subcellular levels. The rapid development of sequencing-based platforms (e.g., [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0), [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq), [Pixel-seq](https://www.cell.com/cell/fulltext/S0092-8674(22)01367-8), [10x Visium HD](https://www.10xgenomics.com/platforms/visium)) and imaging-based platforms (e.g., [10x Xenium](https://www.10xgenomics.com/platforms/xenium), [Vizgen MERSCOPE](https://vizgen.com/merscope-ultra/), [CosMx SMI](https://nanostring.com/products/cosmx-spatial-molecular-imager)) has led to an explosion in the production of ST datasets across diverse tissues and species. 

However, the rapid diversification of spatial transcriptomics platforms has led to heterogeneous data formats, inconsistent metadata schemas, and platform-specific analysis pipelines, hindering standardized ingestion, cross-platform integration, and reproducible analysis. Moreover, while modern ST technologies generate rich, multi-modal spatial data at cellular and subcellular resolution, existing tools lack the scalability and interoperability required for integrated visualization and analysis across modalities and samples—highlighting a critical need for a unified, scalable ecosystem.

To address this challenge, we introduce a scalable, harmonized, and cloud-friendly ecosystem for spatial omics data across platforms, composed of two components:

- **CartLoader** – a harmonized, scalable pipeline for high‑resolution spatial omics data processing and analysis, while retaining the original resolution.
- [**CartoStore**](https://github.com/seqscope/cartostore) – an open-access, cloud-hosted repository for dataset sharing and visualization.
- [**CartoScope**](placeholder): Provides interactive, pixel-level exploration and analysis of spatial omics data.

Together, this system should provide a unified solution for working with raw pixel-level spatial omics data—without sacrificing resolution or reproducibility.

---

## `CartLoader`: A Scalable Spatial Transcriptomics Pipeline

This document introduces `CartLoader`, the core engine of our ecosystem. It provides a modular, reproducible tool to harmonize, integrate, analyze, and visualize raw high-resolution ST data across platforms. Usage examples are provided in the [Vignettes](./vignettes/intro.md). A detailed description (action, inputs, outputs, and parameters) is provided for each module in the [Reference Pages](./reference/intro.md).

**Key Features**

- **Cross-Platform Format Conversion**: Converts raw SGE data from diverse ST platforms into a unified format, enabling consistent downstream processing.
- **Spatial Factor Inference**: Applies [FICTURE](https://www.nature.com/articles/s41592-024-02415-2) to infer spatial factors directly from pixel-level data, capturing biological patterns without requiring cell segmentation.
- **Multi-Modal Alignment**: Aligns and overlays histology images with SGE data so all layers share a common coordinate system for pixel-accurate comparisons.
- **Cloud-Friendly Outputs**: Produces compact, geospatially-indexed data formats suitable for web visualization and cloud storage.
- **Batch Integration and Sample Stitching**: Supports joint analysis and SGE stitching across samples or platforms to reveal shared or differential features across tissues.
- **Modular and Reproducible Workflow**: Orchestrates all steps through a Makefile-based system to ensure scalability, transparency, and reproducibility.

---

## [CartoStore](https://github.com/seqscope/cartostore): Cloud Access and Visualization

As a natural companion to `CartLoader`, [`CartoStore`](https://github.com/seqscope/cartostore) hosts datasets and outputs from `CartLoader` for public access and exploration.

Designed for scalability and interoperability, [`CartoStore`](https://github.com/seqscope/cartostore) uses spatially indexed formats like PMTiles to support interactive visualization and seamless integration.

Explore the [`CartoStore` documentation](https://github.com/seqscope/cartostore) to learn more about dataset access and how to contribute.

---

## [CartoScope](placeholder): Interactive Exploration at Original Resolution

CartoScope is the front‑end application that enables interactive, pixel‑level exploration of spatial omics data without requiring programming expertise.

### Core Capabilities

- **Multi‑layer interactive maps (spatial factors, genes, histology, boundaries)**
- **Pixel‑level visualization preserving subcellular and extracellular signals**
- **Region‑of‑Interest (ROI) based comparative analysis**
- **Session bookmarking and shareable exploration states**

---

## Citations:

<!-- * CartLoader+CartoStore: [doi-to-be-added](link-to-be-added) -->
* FICTURE: [doi.org/10.1038/s41592-024-02415-2](https://www.nature.com/articles/s41592-024-02415-2)
