# CosMx SMI Reformatting

## Overview

`reformat_cosmx` converts the raw AtoMx / CosMx SMI flat-file export into the inputs the
CartLoader pipeline consumes. Unlike other imaging platforms, CosMx is not ingested through
[`sge_convert`](./sge_convert.md); instead this step reads three CSVs and writes:

* a transcript TSV (`X`, `Y`, `gene`, `count`, `cell_id`, `Z`) for FICTURE,
* a cell-metadata CSV (`cell_id`, `X`, `Y`) — the `xy` role for cell decode, and
* a polygon CSV (`cell_id`, `vertex_x`, `vertex_y`) — the `boundaries` role.

Global-pixel coordinates are converted to microns (`--um-per-px`, default `0.12028`, i.e. 1/8.3).
Cell ids are formed as `<fov>_<cell_ID>`. Transcripts with `cell_ID = 0` are kept but tagged
`UNASSIGNED`; `System*` control probes are dropped.

It is invoked automatically by the [`cosmx_smi` profile of `run_together`](./platforms/cosmx_smi.md);
run it directly only for standalone reformatting.

---
## Usage

```bash
cartloader reformat_cosmx \
    --tx   <slide>_tx_file.csv.gz \
    --meta <slide>_metadata_file.csv.gz \
    --poly <slide>-polygons.csv.gz \
    --out  OUT/PREFIX
```

Writes `OUT/PREFIX.transcripts.tsv.gz`, `OUT/PREFIX.metadata.csv.gz`, and
`OUT/PREFIX.polygons.csv.gz`. `--meta` and `--poly` are optional; omit them to reformat only
the transcript file.

---
## Key options

| Option | Default | Description |
|--------|---------|-------------|
| `--tx` | *(required)* | Transcript CSV (`*_tx_file.csv[.gz]`) |
| `--meta` | — | Cell-metadata CSV (`*_metadata_file.csv[.gz]`) |
| `--poly` | — | Polygon CSV (`*-polygons.csv[.gz]`) |
| `--out` | *(required)* | Output prefix |
| `--um-per-px` | `0.12028` | Microns per pixel scaling factor |
| `--offset-x`, `--offset-y` | `0.0` | Global offsets (microns) added to all coordinates |
| `--tx-col-*`, `--meta-col-*`, `--poly-col-*` | CosMx standard | Column-name overrides per input file |

Column defaults follow the standard CosMx export (`fov`, `cell_ID`/`cellID`, `target`,
`x_global_px`/`y_global_px`, `CenterX_global_px`/`CenterY_global_px`, `z`). Override them if
your export uses different headers.
