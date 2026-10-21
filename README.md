# CartLoader

## Overview

`CartLoader` is a scalable, harmonized, and cloud-friendly pipeline for spatial omics data.

Spatial transcriptomics (ST) technologies have revolutionized our ability to map molecular features with remarkable resolution—down to cellular and subcellular levels.  While the rapid development of sequencing-based platforms (e.g., Seq-Scope, Stereo-seq, Pixel-seq, 10x Visium HD) and imaging-based platforms (e.g., 10x Xenium, Vizgen MERSCOPE, CosMx SMI) has led to an explosion in the production of ST datasets across diverse tissues and species, this platform diversity has also introduced heterogeneous data formats and metadata schemas, which further obstruct efforts toward standardized data ingestion, cross-platform format conversion, and reproducible analytical workflows.

**`CartLoader`** is a harmonized, scalable pipeline for high-resolution spatial omics data processing and analysis with retaining the original resolution.

## Quick Start

`cartloader run_together` takes a sample from raw platform output to deployable CartoScope assets (SGE conversion, FICTURE spatial factors, cell analyses, background images, and PMTiles packaging) in one resumable command. Try it with the prebuilt Docker image and its bundled starter dataset:

```bash
docker run --rm -v $(pwd)/out:/out hyunminkang/cartloader:latest \
  run_together --platform seqscope \
    --in-transcript /app/data/sge/transcripts.unsorted.tsv.gz \
    --out-dir /out/seqscope_hippo --id seqscope_hippo --width 18 --n-factor 6,12
```

The packaged result is written to `out/seqscope_hippo/cartl/seqscope_hippo/`. See [Get Started](https://seqscope.github.io/cartloader/vignettes/quickstart/intro/) for the full walkthrough, including a local installation.

## Installation

Step-by-step instructions are in the [Installation Guide](./docs/installation.md) ([online](https://seqscope.github.io/cartloader/installation/)).

## Documentation

The full documentation, including guides for each platform and a reference for every command, is at **https://seqscope.github.io/cartloader/**.
