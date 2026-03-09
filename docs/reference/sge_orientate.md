# SGE Reorientation (`sge_orientate`)

## Overview

`sge_orientate` rotates and flips an SGE to match the desired image orientation.

---
## Requirements

- Input SGE in the unified format (from [SGE format conversion](./sge_convert.md))

---
## Actions

Rotations are clockwise; flips are applied after rotation.

---
## Example Usage

Below is an example with 90 degree rotation and a horizontal flip.

```bash
OUT_PREFIX=r90_hflip                # Replace r90_hflip with your prefix
cartloader sge_orientate  \
    --in-transcript /path/to/transcripts.unsorted.tsv.gz \
    --in-feature  /path/to/feature.clean.tsv.gz \
    --in-minmax /path/to/coordinate_minmax.tsv \
    --out-transcript /path/to/${OUT_PREFIX}.transcripts.unsorted.tsv.gz \
    --out-feature /path/to/${OUT_PREFIX}.feature.clean.tsv.gz \
    --out-minmax  /path/to/${OUT_PREFIX}.coordinate_minmax.tsv \
    --rotate 90 \
    --flip-horizontal
```

---
## Parameters

### Input/Output Parameters
- `--in-transcript` (str): Path to input transcript SGE TSV (gzipped).
- `--in-feature` (str): Path to input feature TSV (gzipped).
- `--in-minmax` (str): Path to input coordinate min/max TSV.
- `--out-transcript` (str): Output path for oriented transcript SGE TSV (gzipped).
- `--out-feature` (str): Output path for feature TSV (copied from input).
- `--out-minmax` (str): Output path for oriented min/max TSV.

### Key Parameters
- `--rotate` (90|180|270): Rotate clockwise (applied before flips).
- `--flip-vertical` (flag): Flip vertically (flipped along the Y-axis; applied after rotation).
- `--flip-horizontal` (flag): Flip horizontally (flipped along the X-axis; applied after rotation).
- `--chunk-size` (int, default: 1_000_000): Rows per chunk when processing the transcript TSV.

---
## Output

An output SGE matrix after orientation.
