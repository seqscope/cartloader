# SGE Reorientation (`sge_orientate`)

## Overview

`sge_orientate` rotates and flips an SGE to match the desired image orientation.

!!! tip "Not part of `run_together`"
    [`run_together`](./run_together.md) does not reorient the SGE; [`sge_stitch`](./sge_stitch.md) calls `sge_orientate` for tiles that need it. Use this page to reorient an SGE directly in step-by-step or custom workflows.

---
## Requirements

- Input SGE in the unified format (from [SGE format conversion](./sge_convert.md)). The coordinate columns must be named `X` and `Y` (fixed; not configurable), and the min/max file must be the two-column `xmin`/`xmax`/`ymin`/`ymax` TSV.

---
## Actions

Rotations are clockwise; flips are applied after rotation. The SGE is transformed within its own bounding box, so the coordinate range is preserved (with X and Y ranges swapped for 90°/270°).

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

All six are required.

- `--in-transcript` (str): Input transcript SGE TSV (gzipped).
- `--in-feature` (str): Input feature TSV (gzipped).
- `--in-minmax` (str): Input coordinate min/max TSV.
- `--out-transcript` (str): Output oriented transcript SGE TSV (gzipped).
- `--out-feature` (str): Output feature TSV (an unchanged copy of the input).
- `--out-minmax` (str): Output oriented min/max TSV.

### Key Parameters
- `--rotate` (`90`|`180`|`270`): Rotate clockwise (applied before flips). Omit for no rotation.
- `--flip-vertical` (flag): Flip vertically (along the Y-axis; applied after rotation).
- `--flip-horizontal` (flag): Flip horizontally (along the X-axis; applied after rotation).
- `--chunk-size` (int): Rows per chunk when processing the transcript TSV (default: 1,000,000).

Equivalent combinations are normalized first (e.g. `--rotate 180` is the same as both flips).

---
## Output

The oriented transcript TSV and min/max TSV, plus a copy of the feature TSV, at the three `--out-*` paths.
