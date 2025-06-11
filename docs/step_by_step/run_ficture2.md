# Factor Inference Analysis

## Overview
After harmonization, Cartloader applies [**FICTURE (Factor Inference of Cartographic Transcriptome at Ultra-high REsolution)**](https://www.nature.com/articles/s41592-024-02415-2) (Si et al., *Nature Methods*, 2024), a segmentation-free spatial factor analysis method that models gene expression as a continuous tissue field. 

By default, FICTURE applies a standard latent Dirichlet allocation (LDA) model on the hexgonal grids of the spatial data to identify **spatial factors** that reflect underlying transcriptional programs, such as cell types, cell states, extracellular structures or gradients. These factors are then projected at the original spatial resolution, retaining the high-fidelity mapping of biological signal across tissue. Alternatively, the factors can be obtained from external single-cell or single-nucleus RNA-seq references datasets, or from other spatially agnostic methods such as Seurat18 and Scanpy27 that learn factors based on grids. FICTURE can use these external factors as priors to initialize the LDA factors or use them directly as input factors (Fig. 1c).



When available, Cartloader also aligns **histology images** to the SGE data. This co-registration results in a unified spatial coordinate system across:

- Raw gene expression
- Computed spatial factors (from FICTURE)
- Tissue morphology (from histology)

This integration allows **pixel-level comparison** across molecular and morphological layers, enabling new insights into spatial organization and gene activity.
