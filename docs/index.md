# A scalable, harmonized, and cloud-friendly ecosystem for spatial transcriptomics data

## Overview
Spatial transcriptomics (ST) technologies have revolutionized our ability to map molecular features with remarkable resolution—down to cellular and subcellular levels. 

While the rapid development of sequencing-based platforms (e.g., [Seq-Scope](https://www.nature.com/articles/s41596-024-01065-0), [Stereo-seq](https://www.bgi.com/global/service/spatial-transcriptome-stereo-seq), [Pixel-seq](https://www.cell.com/cell/fulltext/S0092-8674(22)01367-8), [10x Visium HD](https://www.10xgenomics.com/platforms/visium)) and imaging-based platforms (e.g., [10x Xenium](https://www.10xgenomics.com/platforms/xenium), [Vizgen MERSCOPE](https://vizgen.com/merscope-ultra/), [CosMx SMI](https://nanostring.com/products/cosmx-spatial-molecular-imager)) has led to an explosion in the production of ST datasets across diverse tissues and species, this platform diversity has also introduced heterogeneous data formats and metadata schemas, which further obstruct efforts toward standardized data ingestion, cross-platform harmonization, and reproducible analytical workflows.

To address this challenge, we introduce a scalable, harmonized, and cloud-friendly ecosystem for spatial transcriptomics data across platforms, composed of two components:

- **Cartloader** – a harmonized, scalable pipeline for high-resolution spatial omics data processing and analysis with retaining the original resolution.
- **Cartostore** – an open-access, cloud-hosted repository for dataset sharing and visualization.

Together, this system should provide a unified solution for working with raw pixel-level spatial omics data—without sacrificing resolution or reproducibility.

---

## Cartloader: A Scalable Spatial Transcriptomics Pipeline

This document introduces `cartloader`, which is the core engine of our ecosystem. It provides a modular, reproducible tool to harmonize, integrate, analyze, and visualize raw high-resolution ST data across platforms. The detailed descriptions of its workflow, required inputs, generated outputs, and illustrative examples are provided in the sections of [Step-by-Step](./step_by_step/intro.md) and [Examples]().

### Key Features

- **Cross-Platform Harmonization**: Converts raw spatial gene expression (SGE) data from diverse ST platforms into a unified format, enabling consistent downstream processing.
- **Spatial Factor Inference**: Applies [FICTURE](https://www.nature.com/articles/s41592-024-02415-2) to infer spatial factors directly from pixel-level data, capturing biological patterns without requiring cell segmentation.
- **Multi-Modal Alignment**: Aligns and overlays histology images to SGE data so all layers (histology, SGE, and histological images) share a common coordinate system for pixel-accurate comparisons.
- **Cloud-Friendly Outputs**: Produces compact, geospatially-indexed data formats suitable for web visualization and cloud storage.
- **Batch Integration and Sample Stitching**: Supports joint analysis and SGE stitching across samples or platforms to reveal shared or differential features across tissues.
- **Modular and Reproducible Workflow**: Orchestrates all steps through a Makefile-based system to ensure scalability, transparency, and reproducibility.

---

## Cartostore: Cloud Access and Visualization

As a natural companion to `cartloader`, **cartostore** hosts the output of processed datasets for public access and exploration. Designed for scalability and interoperability, cartostore uses spatially indexed formats like **PMTiles** to support interactive visualization and seamless integration.

Explore the [Cartostore documentation](git@github.com:seqscope/cartostore.git) to learn more about dataset access and how to contribute.

---

## Citations:

<!-- * cartloader+cartostore: [doi-to-be-added](link-to-be-added) -->
* FICTURE: [doi.org/10.1038/s41592-024-02415-2](https://www.nature.com/articles/s41592-024-02415-2)

