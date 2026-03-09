# What is `FICTURE` and `punkst`?

## FICTURE 
[`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2) is a segmentation-free method. It reconstructs the fine-scale tissue structure by first decomposing gene expression patterns across the tissue section into spatial factors and then assigns each pixel to these factors using local context. Biologically, these inferred factors may correspond to specific cell types, functional or physiological states, subcellular domains, or extracellular transcriptomic signatures.

By default, [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2) learns the spatial factors by implementing a standard latent Dirichlet allocation (LDA) model on a hexagonal grid overlay of the spatial coordinates. Optionally, spatial factors can also be derived from external sources, such as single-cell or single-nucleus RNA-seq reference datasets, or from spatially agnostic factor learning methods (e.g., Seurat, Scanpy).


## Punkst

To efficiently run FICTURE-based inference, `CartLoader` integrates [`punkst`](https://github.com/Yichen-Si/punkst), an optimized implementation of FICTURE that maintains output equivalence while enhancing computational scalability and performance. 

Currently, `run_ficture2` uses the `punkst` version of `FICTURE`.
