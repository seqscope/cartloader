# Reference Overview

This section documents the `CartLoader` CLI. Start with [`run_together`](./run_together.md), the recommended end-to-end workflow; the individual commands it runs are documented below for step-by-step or custom use.

## run_together

- [Overview](./run_together.md): `run_together` runs ingest → FICTURE → cell analyses → packaging → images → (optional) annotation and upload, for one sample or many.
- [Specifying Inputs](./run_together_inputs.md): Single-sample CLI, multi-sample sheet, and full JSON config.
- [Image Modalities](./run_together_images.md): Declaring morphology and histology images to package.
- [Platforms (Overview & Status)](./supported_platforms.md): Built-in platform profiles and their status.
    - [10x Xenium](./platforms/xenium.md), [10x Visium HD](./platforms/visium_hd.md), [CosMx SMI](./platforms/cosmx_smi.md), [MERSCOPE](./platforms/merscope.md), [Illumina StrataMap](./platforms/illumina.md), [BGI Stereo-seq](./platforms/stereoseq.md), [Seq-Scope](./platforms/seqscope.md), [Generic](./platforms/generic.md): Expected inputs and example commands per platform.

## Individual Commands

### SGE Preparation

- [SGE Format Conversion](./sge_convert.md): `sge_convert` standardizes raw platform outputs into the unified SGE format, with optional visualization and density filtering.
- [CosMx Reformatting](./reformat_cosmx.md): `reformat_cosmx` turns a CosMx SMI flat-file export into transcript, cell, and boundary inputs.
- [SGE Stitching](./sge_stitch.md): `sge_stitch` stitches two or more SGE tiles into one.
- [SGE Reorientation](./sge_orientate.md): `sge_orientate` rotates and flips an SGE.

### Analysis

- [FICTURE Analysis](./run_ficture2.md): `run_ficture2` trains spatial factors and decodes pixels for one sample.
- [Multi-Sample FICTURE Analysis](./run_ficture2_multi.md): `run_ficture2_multi` trains joint factors across samples and decodes each sample.
- [Multi-Sample Cell Analysis](./run_ficture2_multi_cells.md): `run_ficture2_multi_cells` decodes cell analyses with the joint model.
- [Feature Customization](./feature_customization.md): Filtering features by list, substring, regex, or type at each stage.

### Packaging

- [Asset Packaging](./run_cartload2.md): `run_cartload2` packages SGE and FICTURE outputs into PMTiles and a `catalog.yaml`.
- [Multi-Sample Asset Packaging](./run_cartload2_multi.md): `run_cartload2_multi` packages each sample of a joint run plus a shared `multi-catalog.yaml`.
- [Background Image Import](./import_image.md): `import_image` converts an image into a georeferenced PMTiles layer and registers it in the catalog.
- [Cell Analysis Import](./import_cell.md): `import_xenium_cell` and `import_visiumhd_cell` import platform cell analyses as PMTiles layers.
- [Square Analysis Import](./import_square.md): `import_visiumhd_square` imports Visium HD square-bin analyses as PMTiles layers.

### Upload

- [AWS Upload](./upload_aws.md): `upload_aws` uploads catalog-listed files to an S3 path.
- [Zenodo Upload](./upload_zenodo.md): `upload_zenodo` uploads to a new or existing Zenodo deposition.

### AI Annotation

- [Annotate CartLoad Folder](./anno_cartload_folder.md): `anno_cartload_folder` AI-annotates every factor in a catalog and records the alias files.
- [Annotate Factors with an LLM](./annotate_factors_with_llm.md): `annotate_factors_with_llm` labels all factors of a model in one LLM call and writes an interactive report.

## Legacy Orchestrators

- [Overview](./legacy.md): Status of the legacy platform orchestrators, what `run_together` replaces, and known issues.
- [run_xenium](./run_xenium.md): `run_xenium`, the legacy single-sample Xenium orchestrator.
- [Xenium Add-ons](./xenium_addon.md): `load_xenium_ranger` writes a JSON manifest of Xenium Ranger outputs.
- [run_visiumhd](./run_visiumhd.md): `run_visiumhd`, the legacy single-sample Visium HD orchestrator.
- [Visium HD Add-ons](./visiumhd_addon.md): `load_space_ranger` writes a JSON manifest of Space Ranger outputs.
