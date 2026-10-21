# CosMx SMI Reformatting

## Overview

`reformat_cosmx` converts the raw AtoMx / CosMx SMI flat-file export into the inputs the
CartLoader pipeline consumes. Unlike other imaging platforms, CosMx is not ingested through
[`sge_convert`](./sge_convert.md) by `run_together`; instead this step reads three CSVs and writes:

* a transcript TSV (`X`, `Y`, `gene`, `count`, `cell_id`, `Z`, plus `CellComp` when the input has it) for FICTURE,
* a cell-metadata CSV (`cell_id`, `X`, `Y`) — the `xy` role for cell decode, and
* a polygon CSV (`cell_id`, `vertex_x`, `vertex_y`) — the `boundaries` role.

Global-pixel coordinates are converted to microns (`--um-per-px`, default `0.12028`, i.e. 1/8.3).
Cell ids are formed as `<fov>_<cell_ID>`. Transcripts with `cell_ID = 0` are kept but tagged
`UNASSIGNED`; `System*` control probes are dropped. Spaces in gene names are replaced with `-`.
If the transcript CSV has a `CellComp` column (the subcellular compartment, e.g. `Nuclear`, `None`; present
in some exports only), it is kept as a trailing `CellComp` column; rename the input column with `--tx-col-cellcomp`.

!!! tip "Usually run for you"
    The [`cosmx_smi` profile of `run_together`](./platforms/cosmx_smi.md) calls `reformat_cosmx` as its ingest stage. Use this page to run it directly in step-by-step or custom workflows.

---
## Usage

```bash
cartloader reformat_cosmx \
    --tx   <slide>_tx_file.csv.gz \
    --meta <slide>_metadata_file.csv.gz \
    --poly <slide>-polygons.csv.gz \
    --out  OUT/PREFIX \
    --exclude-feature-regex '^(Negative|NegPrb)'
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
| `--out-suffix-tx` | `.transcripts.tsv.gz` | Suffix of the output transcript TSV |
| `--out-suffix-meta` | `.metadata.csv.gz` | Suffix of the output cell-metadata CSV |
| `--out-suffix-poly` | `.polygons.csv.gz` | Suffix of the output polygon CSV |
| `--um-per-px` | `0.12028` | Microns per pixel scaling factor |
| `--offset-x`, `--offset-y` | `0.0` | Global offsets (microns) added to all coordinates |
| `--tx-col-*`, `--meta-col-*`, `--poly-col-*` | CosMx standard | Column-name overrides per input file |

Column defaults follow the standard CosMx export (`fov`, `cell_ID`/`cellID`, `target`,
`x_global_px`/`y_global_px`, `CenterX_global_px`/`CenterY_global_px`, `z`). Override them if
your export uses different headers.

### Feature filtering

Applied while the transcript TSV is written, so a dropped feature is absent from every downstream product of this ingest. The metadata and polygon files are not affected.

| Option | Description |
|--------|-------------|
| `--include-feature-list` | File listing the feature names to keep (one per line) |
| `--exclude-feature-list` | File listing the feature names to drop (one per line) |
| `--include-feature-regex` | Keep only features whose name matches |
| `--exclude-feature-regex` | Drop features whose name matches |

A feature is kept only if it passes every filter given; an exclusion wins over an inclusion. Regexes use **substring (search)** semantics (anchor with `^`/`$` as needed) and are matched against the gene name after spaces are replaced with `-`. See [Feature Customization](./feature_customization.md) for how filters behave across commands.
