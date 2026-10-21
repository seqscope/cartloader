# SGE Stitching (`sge_stitch`)

## Overview

`sge_stitch` stitches two or more SGE tiles into a single SGE, with optional per-tile orientation.
It can also run density-based filtering on the stitched SGE.

!!! tip "Not part of `run_together`"
    [`run_together`](./run_together.md) does not stitch tiles. Use `sge_stitch` in step-by-step or custom workflows; its stitched transcript TSV can then be given to `run_together` as an already-ingested `transcript` (`--in-transcript`; see [Transcript sources](./run_together_inputs.md#transcript-sources)).

---
## Requirements

- Input SGE tiles in the unified format (from [SGE format conversion](./sge_convert.md)), with `X` and `Y` coordinate columns.
- `cartloader` on the `PATH` (the generated Makefile calls `cartloader` sub-commands).
- `gzip`; the `ficture` Python package (installed with CartLoader) for `--filter-by-density`.

---
## Example Usage

Define each input tile as a comma-separated tuple of paths and layout/orientation fields:

```bash
cartloader sge_stitch \
    --in-tiles \
        "/path/to/tileA/transcripts.unsorted.tsv.gz,/path/to/tileA/feature.clean.tsv.gz,/path/to/tileA/coordinate_minmax.tsv,0,0,,false,false" \
        "/path/to/tileB/transcripts.unsorted.tsv.gz,/path/to/tileB/feature.clean.tsv.gz,/path/to/tileB/coordinate_minmax.tsv,0,1,90,false,true" \
    --out-dir /path/to/output/dir \
    --units-per-um 1.0 \
    --precision 2
```

Tile A is placed as is (empty `rotate` field); tile B is rotated 90° clockwise, then flipped horizontally.

---
## Parameters

### SGE Stitch
- `--in-tiles` (list of str): One tuple per input tile (see below).
- `--out-dir` (str, required): Output directory for the stitched (and filtered) SGE and the Makefile.

!!! info "Input Tile Tuple Format"

    Each tuple has exactly eight comma-separated fields; quote each tuple:
    ```
    "transcript_path,feature_path,minmax_path,row,col,rotate,vertical_flip,horizontal_flip"
    ```

    - `transcript_path` (str): The tile's transcript-indexed SGE TSV (`.tsv.gz`).
    - `feature_path` (str, optional): The tile's per-gene count TSV. Leave empty (or `None`) to generate it with `cartloader sge_supp`.
    - `minmax_path` (str, optional): The tile's coordinate min/max TSV. Leave empty (or `None`) to generate it with `cartloader sge_supp`.
    - `row`, `col` (int): Grid position of the tile. Each grid cell is as wide as the widest tile in its column and as tall as the tallest tile in its row, and each tile is centered in its cell; `col` increases along X and `row` along Y. A `(row, col)` pair may be used only once.
    - `rotate`: Clockwise rotation, `90`, `180`, or `270`; leave the field **empty** for no rotation (`0` is not accepted).
    - `vertical_flip`, `horizontal_flip`: `true`/`false` (case-insensitive); flips are applied after rotation.

    Generated feature/min-max files are written **next to the tile's transcript file**: `<id>.feature.tsv.gz` and `<id>.minmax.tsv`, where `<id>` is the transcript file name without `.tsv.gz` and `.transcripts`; or `feature.clean.tsv.gz` and `coordinate_minmax.tsv` when the transcript is named `transcripts.unsorted.tsv.gz`.

!!! warning "Known issues"
    - A tile whose orientation reduces to flips only — a flip without rotation, or `180` with or without a single flip — currently fails at its `sge_orientate` step. Pre-orient such tiles with [`sge_orientate`](./sge_orientate.md) and pass the result with an empty `rotate` field and both flips `false`.
    - `--sge-visual` currently fails while the Makefile is being generated; omit it.

??? note "Auxiliary SGE Stitching Parameters"
    We recommend using default values; override only if needed.

    **Auxiliary Output Parameters**:

    - `--out-transcript` (str): Output transcript-indexed SGE under `--out-dir` (default: `transcripts.unsorted.tsv.gz`).
    - `--out-feature` (str): Output per-gene count TSV under `--out-dir` (default: `feature.clean.tsv.gz`).
    - `--out-minmax` (str): Output coordinate min/max TSV under `--out-dir` (default: `coordinate_minmax.tsv`).
    - `--out-tile-minmax` (str): Output per-tile min/max and offset TSV under `--out-dir` (default: `coordinate_minmax_per_tile.tsv`).
    - `--out-json` (str): Output JSON manifest of the SGE (default: `<out-dir>/sge_assets.json`).

    **Auxiliary Colname Parameters**:

    - `--colname-count` (str): Count column (default: `count`).
    - `--colname-feature-name` (str): Feature name column (default: `gene`).
    - `--colname-feature-id` (str): Feature ID column, if present (default: none).
    - `--colname-x`, `--colname-y` (str): X/Y columns (defaults: `X`, `Y`).
    - `--colnames-others` (list): Additional columns to keep.

    **Auxiliary Scale and Precision Parameters**:

    - `--units-per-um` (float): Coordinate units per µm of the input tiles (default: 1.0). Tiles from `sge_convert` are already in µm; keep the default.
    - `--precision` (int): Decimal places of the stitched coordinates (default: 2).

### Density-based Filtering

- `--filter-by-density` (flag): Enable density-based filtering of the stitched SGE.
- `--out-filtered-prefix` (str): Prefix of the filtered outputs (default: `filtered`).

??? note "Auxiliary Density-based Filtering Parameters"
    We recommend using default values; override only if needed.

    - `--radius` (int): Radius for the polygon area calculation (default: 15).
    - `--quartile` (int): Quartile for the polygon area calculation (default: 2).
    - `--hex-n-move` (int): Sliding step (default: 1).
    - `--polygon-min-size` (int): Minimum polygon size (default: 500).

### Run Parameters

- `--dry-run` (flag): Write the Makefile and print the commands without executing.
- `--restart` (flag): Ignore existing outputs and rerun all steps.
- `--n-jobs`, `-j` (int): Number of parallel jobs (default: 1).
- `--makefn` (str): Makefile name under `--out-dir` (default: `sge_stitch.mk`).

### Environment Parameters

- `--gzip` (str): Path to `gzip` (default: `gzip`).

---
## Output

Under `--out-dir`:

- The stitched SGE: `transcripts.unsorted.tsv.gz`, `feature.clean.tsv.gz`, `coordinate_minmax.tsv` (see the [unified SGE format](./sge_convert.md#unified-sge-matrix)).
- `coordinate_minmax_per_tile.tsv`: per tile, the local min/max, the X/Y offsets applied, and the resulting global bounds in µm.
- `feature.distribution.tsv.gz`: per feature, the number of tiles containing it and its count in each tile.
- `sge_assets.json`: paths of the SGE to use downstream, in the same format as [`sge_convert`'s](./sge_convert.md#sge-manifest-sge_assetsjson); written after the Makefile runs.
- Oriented tiles (for tiles with a rotation): `r<row>c<col>.rot<deg>[.vfip][.hfip].transcripts.unsorted.tsv.gz`, with matching `.feature.clean.tsv.gz` and `.coordinate_minmax.tsv`.
- With `--filter-by-density`: `<prefix>.transcripts.unsorted.tsv.gz`, `<prefix>.feature.lenient.tsv.gz`, `<prefix>.coordinate_minmax.tsv` (and boundary GeoJSONs); `sge_assets.json` then points to these files.
- `sge_stitch.done`, `sge_density_filtering.done`, and per-tile `.done` flags; the Makefile `sge_stitch.mk`.
