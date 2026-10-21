# Specifying Inputs to `run_together`

`run_together` accepts inputs in **three ways**, from most convenient to most expressive. They compose — you can set the base run on the CLI and add only the extras in a config JSON.

| Mode | Flag | Use when |
|------|------|----------|
| **1. Single sample (direct CLI)** | `--in-dir` / `--in-prefix` **or** explicit `--in-*` file flags | One sample; a standard platform directory, or a few individually-named files. |
| **2. Multi-sample (sample sheet)** | `--samples samples.tsv` | Several samples (a joint model by default). A wide TSV, one row per sample. |
| **3. Full config (JSON)** | `--config run.json` | Per-sample settings differ, or you need analyses/images beyond the profile defaults. |

Samples from all three are **combined**: the config's `samples`, then the `--samples` rows, then the one sample built from the single-sample flags. All samples sharing one `--out-dir` form one joint run.

Whichever you use, the **platform profile** (`--platform <name>`) supplies the defaults: which files to auto-detect, the default FICTURE and cell analyses, and the image conventions. See [Supported Platforms](./supported_platforms.md) for what each profile expects and the per-platform pages for worked examples.

---
## Mode 1 — Single sample (direct CLI)

Two shapes, depending on how your files are laid out.

**A standard platform directory** — point `--in-dir` at it and the profile auto-detects every role inside:

```bash
cartloader run_together --platform 10x_xenium \
    --in-dir /data/sample1 --out-dir OUT \
    --width 12 --n-factor 24
```

**Individually-named files** (not in a standard directory) — name each explicitly. Each flag is the CLI equivalent of a sample-sheet column:

```bash
cartloader run_together --platform merfish \
    --raw-transcript   /data/molecules.csv \
    --in-cell-xy       /data/cell_metadata.csv \
    --in-cell-boundary /data/cell_boundaries.csv \
    --id MYSAMPLE --out-dir OUT
```

| Flag | Sample-sheet column | Meaning |
|------|--------------------|---------|
| `--in-dir` | `in_dir` | Raw platform directory → roles auto-detected inside it. A platform may instead consume the directory whole (Seq-Scope: it *is* the MEX directory). |
| `--in-prefix` | `in_prefix` | Raw platform **path prefix** → roles auto-detected by suffix (Stereo-seq). |
| `--raw-transcript` | `raw_transcript` | Raw transcript CSV/TSV (or `.parquet`, converted first) to ingest through `sge_convert`. |
| `--in-transcript` | `transcript` | An **already-ingested** transcript TSV; skips ingest. Mutually exclusive with `--raw-transcript`. |
| `--in-cell-xy` | `xy` | Cell centroids / metadata file. |
| `--in-cell-boundary` | `boundaries` | Cell boundary polygons. |
| `--in-cellxgene` | `cellxgene` | Cell×gene matrix CSV → converted to a MEX (drives cell clustering). |
| `--id` | `id` | Sample id. With `--out-dir` it defaults to **`rep1`** (the packaged directory and catalog id become `<out-dir basename>-rep1`); with `--out-root`, to the `in_dir` / `in_prefix` basename. |

`--in-dir` is just a one-row sample sheet with only `in_dir`; the explicit `--in-*` flags are a one-row sheet with those columns. Other roles (`clusters`, `mex`, `cell_tsv`, `gef`, …) have no single-sample flag: use a one-row sample sheet or a config (or a ready-made [`--list-*` file](#supplying-your-own-list-files)). Images use [`--image`](./run_together_images.md#inline-cli-image).

---
## Mode 2 — Multi-sample (sample sheet)

A joint multi-sample run is a single-sample run with `--samples samples.tsv` in place of `--in-dir`. All rows share one `--out-dir`, so they train **one joint FICTURE model**. For the common case you only need `id` + `in_dir`:

```bash
cartloader run_together --platform 10x_xenium \
    --samples samples.tsv --out-dir OUT --width 12 --n-factor 24 -j 8
```

```
id    in_dir
s1    /data/s1
s2    /data/s2
s3    /data/s3
```

For **independent** per-sample models, use `--out-root` instead of `--out-dir`; each row becomes its own `<out_root>/<id>/`. A sheet with a single row and `--out-dir` is a single-sample run.

### The sample sheet is a wide table of input roles

Columns map to per-sample **input roles**. Every column is optional except that each sample needs a [transcript source](#transcript-sources). Column names are case-sensitive; the listed **aliases** are accepted interchangeably.

| Column (aliases) | Role | Meaning |
|---|---|---|
| `id` | — | Sample identifier. Defaults to the `in_dir` / `in_prefix` basename; a lone `--out-dir` sample without an `id` is always `rep1`. |
| `in_dir` | — | Raw platform directory → the profile **auto-detects** roles inside it. |
| `in_prefix` | — | Raw platform path prefix → the profile auto-detects roles by **suffix** (`{prefix}.tissue.gef`, …); see [BGI Stereo-seq](./platforms/stereoseq.md). |
| `raw_transcript` | — | A **raw** transcript file to ingest (e.g. MERSCOPE `detected_transcripts.csv`, Xenium `transcripts.csv.gz`) through `sge_convert`. |
| `transcript` (`tsv`) | transcript | A **pre-converted** transcript TSV (the *output* of ingest) → **skips ingest**. Not for raw CSVs. |
| `xy` (`cell_xy`) | xy | Cell centroids file. |
| `boundaries` (`cell_boundary`, `cell_boundaries`) | boundaries | Cell boundary polygons. |
| `clusters` | clusters | External cluster labels. |
| `mex` (`mex_dir`) | mex | MEX directory. Or give the explicit triple `mex_bcd` / `mex_ftr` / `mex_mtx` when one directory does not apply. |
| `cellxgene` (`cell_by_gene`) | cellxgene | Cell×gene matrix CSV (e.g. MERSCOPE `cell_by_gene.csv`) → converted to a MEX that drives cell clustering (works without boundaries). Converted only for a sample that goes through ingest, and not on `cosmx_smi` / `stereoseq`. |
| `gef` / `cellbin_gef` | gef, cellbin_gef | Stereo-seq binary GEFs, when they do not match `in_prefix` + the standard suffix. |
| `cell_tsv` | cell_tsv | A standalone pixel TSV with a header row `X  Y  gene  count  cell_id` (columns read by name) that supplies cell counts on its own (Stereo-seq cell bins). |
| `dapi` | — | One DAPI **OME-TIFF** → a colorized `dapi` layer. A name ending in `morphology.ome.tif` (Xenium's z-stack, any prefix) is imported with `--use-middle-page --high-memory`. A value not ending in `.ome.tif`/`.ome.tiff` is imported without OME→PNG conversion and has no georeferencing source; give such images through [`--images`](./run_together_images.md) with `um_per_pixel` instead. |
| `hne` | — | H&E image. On Visium HD it gets the profile's scale handling (see [Visium HD](./platforms/visium_hd.md)); elsewhere it is a plain `rgb` layer, which needs an `.ome.tif` (scale read from its metadata) or, as on Seq-Scope, a profile `image_defaults` that supplies the bounds. Otherwise use `--images` with the georeferencing fields. Used as given (not resolved against `in_dir`) on Visium HD. |
| `tissue` | — | The sample's tissue, used by `--anno-deep` (falls back to `--tissue`). |

Other images are **not** sample-sheet columns — see [Image Modalities](./run_together_images.md). A config sample may also carry its own `images` list (JSON only).

**Rules**

- An explicit column **overrides** auto-detection for that role. The cell values `` (empty), `-`, `.`, and `NA` all mean *unset*.
- Role columns, `dapi`, and (except on Visium HD) `hne` are resolved against `in_dir` when relative and not found as given; `raw_transcript` is not — give it as an absolute path (or relative to the working directory).
- Sheet columns are **not validated** on their own: a misspelled column (e.g. `transcripts`) is ignored. The transcript-source check below lists unrecognized columns when it fails. (Config JSON samples are validated strictly; see [Validation](#validation-of-config-keys).)

**Parquet inputs.** `raw_transcript`, `xy`, `boundaries` and `clusters` may point at a `.parquet` file (Xenium ships `transcripts.parquet`, `cells.parquet`, `cell_boundaries.parquet` beside the `.csv.gz` files). Their consumers read CSV only, so each is converted to `.csv.gz` under `OUT/tsv/<id>/parquet2csv/` with `cartloader parquet_to_csv_rapid` as the first ingest step — announced with a `NOTE:` at planning time — and every later stage uses the converted file. A failed or empty conversion aborts the run.

### Transcript sources

Each sample gets its transcripts from exactly one source:

| Source | Applies to |
|--------|-----------|
| **`transcript`** (`--in-transcript`) — an already-ingested TSV; ingest is skipped | every platform |
| **`raw_transcript`** (`--raw-transcript`) — a raw CSV/TSV ingested by `sge_convert` | platforms that ingest a transcript CSV: `10x_xenium`, `merfish`, `seqscope` (raw per-molecule TSV), `generic` |
| **`in_dir`** — the platform's raw directory | `10x_xenium`, `10x_visium_hd`, `cosmx_smi`, `merfish`, `seqscope` (the MEX directory); `generic` (a `transcripts.tsv.gz` inside it) |
| **`in_prefix`** or a **`gef`** column | `stereoseq` |
| **`mex`** — the MEX directory | `illumina` |

The source is checked **before anything is planned**: a sample with none, or one whose file is missing (including an `in_dir` without the platform's expected input), stops the run with an error naming each affected sample and listing any unrecognized sheet columns. The check is skipped when the ingest stage is not run (`--only` / `--skip ingest`), since a resumed run reuses the already-converted transcripts.

```
# in_dir (auto-detect), a pre-converted TSV, and a raw CSV named explicitly:
id    in_dir            transcript                    raw_transcript                 boundaries
s1    /data/s1
s2                      /data/s2/transcripts.tsv.gz
s3                                                    /data/s3/detected.csv.gz       /data/s3/bounds.csv.gz
```

An already-ingested transcript uses the layout ingest writes (`transcripts.unsorted.tsv.gz`): tab-separated, a header row `X  Y  gene  count` (optionally a fifth `cell_id` column), coordinates in µm.

### Column-name overrides

When input **columns** are non-standard, name them (CLI, or the `ingest` / `roles` blocks of the profile/`--config`). The CLI flags apply to every sample in the run:

- `--colname-transcript-x/-y/-feature/-count` — the raw transcript's coordinate/gene/count columns (`ingest.csv_colnames`). Without a count column each row counts as 1.
- `--colname-xy-cell/-x/-y` — the xy file's columns (`roles.xy`; use `--colname-xy-cell ''` for an unnamed index column).
- `--colname-boundary-cell/-x/-y` — the boundary file's columns as read by `spatula tsv-add-cell-id` (`roles.boundaries`; MERSCOPE, Illumina).

!!! note "Boundary files for polygon rendering"
    The cell-polygon layer is built by `run_cartload2` → `run_boundaries2pmtiles`, which expects a comma-separated file whose first three columns are `cell_id`, `vertex_x`, `vertex_y` (vertices grouped by cell, µm). The `--colname-boundary-*` overrides do not reach this step.

!!! note "Transcript that already carries a `cell_id` column (rare)"
    If your raw transcript CSV already has a per-molecule cell-id column, name it with
    `--colname-transcript-cell <name>` (or `ingest.csv_colname_cell` in `--config`).
    The column is carried through ingest to transcript column 5 and used directly for
    cell analysis, skipping `spatula tsv-add-cell-id`. A `cellxgene` MEX, if also
    present, still drives the clustering.

---
## FICTURE mode (de-novo vs. projection)

Independent of how inputs are specified, the CLI selects the base FICTURE work. Exactly one mode is the base (a config JSON can add more analyses on top).

=== "De-novo training (default)"

    Train new LDA models. `--width` and `--n-factor` accept comma lists → the cross-product is trained (**multiple widths supported**). They override the width / factor count of the profile's training entry; either alone keeps the profile's value for the other.

    ```bash
    cartloader run_together --platform 10x_xenium --in-dir IN --out-dir OUT \
        --width 12,18 --n-factor 24,48
    ```

=== "Projection-only"

    Reuse **already-trained** models — no LDA training runs. Point `--project-models` at one or more **existing FICTURE directories**; `run_together` reads each `ficture.params.json` and re-projects every model it lists onto the current data (at `--width` if given, else each model's training width). The profile's training entry is dropped.

    ```bash
    cartloader run_together --platform 10x_xenium --in-dir IN --out-dir OUT \
        --project-models /prev/run/fic --width 12
    ```

    Ingest still runs (the *data* is current); only training is skipped.

=== "No factor analysis"

    Host a dataset with **no FICTURE at all** — transcripts, the SGE raster, and histology, with no factor layers. `--no-ficture` runs only FICTURE's tiling step (`run_ficture2_multi --prepare-only`), and packaging reads the resulting tiled TSV directly.

    ```bash
    cartloader run_together --platform seqscope --in-dir IN --out-dir OUT --no-ficture
    ```

    Applies to **every** platform. Cell analyses are skipped too (they decode against a model). Cannot be combined with `--n-factor` or `--project-models`. The hexagon files are still built alongside the tiles (at `--width`, else the profile's width, else `12`), so adding factors later re-uses them instead of re-tiling — rerun without `--no-ficture` in the same `--out-dir`.

Cell analyses decode against the model with the **largest factor count** unless an analysis names its own `model_id`.

**Hexagon MEX export** (any mode): `--segment-10x` additionally writes each sample's hexagon file as a 10x MEX directory, `fic/samples/<id>/<id>.hex_<width>.mex/` (`barcodes.tsv.gz` with hexagon centers as `x:y`, `features.tsv.gz`, `matrix.mtx.gz`). It converts the very hexagon files the factor analysis uses (same `--min-ct-per-unit-hexagon` filter), so the two are consistent. By default every analysis' `--width` is exported; `--segment-width-10x 12,24` picks the widths instead (extra widths get their hexagon files built too). Works with `--no-ficture` — the way to get hexagon MEX files without any factor analysis, e.g. for Seq-Scope:

```bash
cartloader run_together --platform seqscope --in-dir IN --out-dir OUT --no-ficture \
    --segment-10x --segment-width-10x 12,24
```

Config equivalent: top-level `"segment_10x": true` or `{ "widths": "12,24" }` (a CLI flag wins). The directories are recorded under `mex` in the per-sample `ficture.params.json` and in `ficture.multi.params.json`; see [`run_ficture2_multi --segment-10x`](./run_ficture2_multi.md#actions).

**Tiling** (run-wide, rarely needed): `--tile-size` and `--tile-buffer` are forwarded to `run_ficture2_multi` (punkst `multisample-prepare`) and otherwise keep its defaults of `500` µm and `1000` lines. punkst rejects a tile size under 20× the hexagon *side* length (width / √3) and recommends 50–100×, so only a wide hexagon needs this — e.g. `--width 100` wants roughly `--tile-size 5000`. The buffer is a per-tile line buffer, not a spatial size. Config equivalents: top-level `"tile_size"` / `"tile_buffer"` (a CLI flag wins).

**Common decode overrides** (else profile / built-in default): `--exclude-feature-regex`, `--include-feature-list` / `--exclude-feature-list`, the `--ingest-*-feature-*` family, `--min-ct-per-unit-hexagon` (default `50`), `--min-ct-per-unit-train`, `--cell-min-cell-count` / `--cell-min-feature-count`, and `--always-single-molecule` / `--never-single-molecule` (default: single-molecule **ON** for pixel FICTURE, **OFF** for cell decode). For these, an explicit CLI flag wins over a `--config`/profile value, which wins over the built-in default. Settings that do not follow that order: `resources.n_jobs` / `resources.threads` in a config win over `-j` / `--threads` (see [Mode 3](#mode-3-full-config-json)), and `--units-per-um` has no effect on `seqscope` or `10x_visium_hd` (see [Coordinate units](#coordinate-units-at-ingest)).

### Feature filtering: two independent layers

Feature filters come in two flavours that do **not** interact. Each takes a regex and/or a plain text file of feature names (one per line), and each has a default that applies to **every** platform:

| | default regex | rationale |
|---|---|---|
| ingest | <code>^(Unassigned&#124;Neg&#124;BLANK&#124;Blank&#124;Intergenic&#124;Deprecated&#124;System&#124;NCS-&#124;NCP-)</code> | technical artifacts — not genes, so they never enter the data |
| analysis | <code>^(Gm[0-9]&#124;MT-&#124;mt-&#124;Rps&#124;Rpl)</code> | real genes that distort a factorization — kept in the data, kept out of the models |

Pass `''` to either flag to disable its default; a profile or `--config` value overrides it.

**1. Analysis (FICTURE) filters** — `--include-feature-list` / `--exclude-feature-list` / `--exclude-feature-regex`. These restrict the **factor model only** (`lda4hex --features`); the data — counts, pseudobulk, DE — keep every gene. `--include-feature-list` is the natural place for a highly-variable-gene set: the model is built on those genes, but all genes are still emitted afterwards.

| stage | affected |
|---|---|
| ingest, tiled transcripts, feature list, packaged tiles | **no** — every gene is kept |
| pixel FICTURE: the LDA/projection **model** (and hence pixel decoding, which reads that model) | yes |
| cell clustering: the LDA is a projection onto the pixel model, so it inherits the same restriction | yes (via the model) |
| cell **counts, pseudobulk, DE** | **no** — excluded/non-HVG genes stay in and reappear here |

A list and the regex **combine**: the regex narrows what the list leaves. The features the pixel model was trained on are written to `fic/multi.selected_features.tsv` and recorded per model as `feature` in `ficture.params.json`.

Because the restriction lives in the model (not the counts), the cell **pseudobulk and DE report every gene**, and the pixel decode — which reads the restricted model — reports only the model's genes. That asymmetry is intentional: the pixel decode *is* the model, the cell pseudobulk is an independent aggregate of the raw counts.

**2. Ingest filters** — `--ingest-include-feature-list` / `--ingest-exclude-feature-list` / `--ingest-include-feature-regex` / `--ingest-exclude-feature-regex`. These drop features from the transcript TSV as it is written, so a filtered feature is gone from **everything** downstream: the feature list, the packaged tiles, the browser's gene list and every analysis. Use them for features that should not be part of the dataset at all (negative-control probes, blanks); use the analysis filters for features that should be visible but not drive the factorization.

The ingest regex always replaces `sge_convert`'s per-platform default (e.g. Xenium's negative-probe pattern), so one pattern governs every platform.

- Applied by whichever step writes the transcript: `sge_convert`, `reformat_cosmx` (CosMx), and the Stereo-seq cell-bin conversion.
- Also applied to the **cell count matrices** (`mex2sptsv` / `pixel2sptsv`), so a MEX-derived cell source (a `cell_by_gene` matrix or an external `mex` role) — which never passes through `sge_convert` — drops the same technical artifacts. This is the *only* feature filter on the cell counts; the FICTURE filters deliberately are not applied there (see above), so the cell pseudobulk keeps every real gene.
- **Ignored, with a warning, for a sample that supplies an already-ingested `transcript`** — that file is used as given. Filter it beforehand, or supply the raw input.
- On the MEX-based platforms (`10x_visium_hd`, `seqscope`, `illumina`) filtering happens inside `spatula convert-sge`, which accepts only **one** include-type and **one** exclude-type filter; a list and a regex of the same polarity is an error there. Resolve them into one list with `cartloader feature_select`.

!!! warning "Count thresholds and the analysis filter"
    `--min-ct-per-unit-hexagon` is applied over **all** genes, while the model is fit over the **restricted** set, so restricting to a small panel leaves hexagons whose surviving counts are low; `--min-ct-per-unit-train`, counted over the features the model uses, is the threshold to lower. The cell counts are not restricted by the analysis filters (only by the ingest filters), so `--cell-min-cell-count` counts over every gene that survived ingest; it applies to cells built from a transcript or `cell_tsv`, not to MEX-derived cells. Ingest filters run before any counting.

---
## Coordinate units at ingest

`--units-per-um <n>` tells ingest how many coordinate units of the **raw input** make one micron — the factor `sge_convert` divides the input X/Y by. `1` means the input is already in microns, `1000` means nanometers, `2` is Stereo-seq's 500 nm bins. Everything downstream of ingest is in microns, so this is the one place the input's unit convention is declared.

Normally you never set it: each platform's ingest preset knows its own convention. Set it when a dataset departs from that convention — notably **Illumina StrataMap**, which ships two barcode formats:

| Barcode | Units | `--units-per-um` |
|---|---|---|
| `SBC:433503:2393851` (older) | nanometers | `1000` |
| `SBC:686.951:4668.15` (current) | microns | `1` |

Ingest **detects which of the two** a StrataMap barcode file uses (a fractional coordinate means microns; otherwise the coordinate magnitude decides) and reports the value it picked, so both formats run correctly with no flag. `--units-per-um` overrides the detection; if the value you give contradicts the file, the run warns and uses yours.

- Applied by `sge_convert`, in both of its routes (MEX via `spatula convert-sge`, CSV via `sge_format_generic`), so it covers every platform whose ingest is `sge_convert`. On `stereoseq` it also scales the cell-bin TSV, keeping the pixel and cell coordinates on one system.
- **Not supported on `cosmx_smi`**, whose ingest (`reformat_cosmx`) writes the transcript TSV itself. Asking for it there is an error.
- **No effect on `seqscope` from the CLI or `ingest.units_per_um`:** the profile's `ingest.extra_flags` already contains `--units-per-um 1000`, which comes later on the `sge_convert` command line and takes precedence. To change it, replace `ingest.extra_flags` in a config (see the [Seq-Scope page](./platforms/seqscope.md#coordinates-and-count-columns)).
- **No effect on `10x_visium_hd`:** `sge_convert` replaces it with the value derived from `microns_per_pixel` in the scale factors JSON (`--scale-json`) that the profile passes.
- **Ignored, with a warning, for a sample that supplies an already-ingested `transcript`** — that file is used as given, in microns.
- Boundary/centroid files are **not** rescaled by this flag: they must already be in microns (see the platform pages).

```bash
# a current (micron) StrataMap barcode file, stated explicitly
cartloader run_together --platform illumina --samples samples.tsv --out-dir OUT --units-per-um 1
```

The config equivalent lives in the `ingest` block (a CLI `--units-per-um` wins over it):

```jsonc
{ "platform": "illumina", "ingest": { "units_per_um": 1 } }
```

---
## Coordinate jitter at ingest

`--jitter-xy <um>` adds a uniform random offset in `[-um, +um]` to each transcript's X and Y — drawn **independently per transcript and per axis** — as the transcript TSV is written. It is off by default (`0`).

Use it on coarse-resolution platforms whose molecules sit on a lattice rather than at measured positions. On 10x Visium HD, for example, every transcript in a 2 µm bin carries that bin's single coordinate, so the transcript cloud is a grid of stacked points; hexagon binning, pixel decoding and the rendered tiles all see the lattice. Jittering by roughly **half the bin pitch** (`--jitter-xy 0.8` for 2 µm bins) spreads each bin's molecules across the area they came from.

Like the ingest feature filters, this rewrites the data itself: the jittered coordinates are what every later stage reads — hexagons, FICTURE, the packaged tiles, the browser.

- Applied by `sge_convert`, in both of its routes: `spatula convert-sge` for MEX input (`10x_visium_hd`, `seqscope`, `illumina`) and `sge_format_generic` for CSV input (`10x_xenium`, `merfish`, `stereoseq`'s bin1 GEM, `generic`). Jitter is in **microns**, applied after the input's units are converted.
- **Not supported on `cosmx_smi`**, whose ingest (`reformat_cosmx`) writes the transcript TSV itself. Asking for it there is an error, not a silently dropped flag.
- On `stereoseq` it covers the bin1 pixel transcript only; the cell-bin TSV (`convert_stereoseq_cellbin`) has no jitter option and stays on the original grid. The run warns when you ask.
- **Ignored, with a warning, for a sample that supplies an already-ingested `transcript`** — that file is used as given.

```bash
cartloader run_together --platform 10x_visium_hd --in-dir IN --out-dir OUT --jitter-xy 0.8
```

The config equivalent lives in the `ingest` block (a CLI `--jitter-xy` wins over it):

```jsonc
{ "platform": "10x_visium_hd", "ingest": { "jitter_xy": 0.8 } }
```

---
## Mode 3 — Full config (JSON)

Escalate to `--config run.json` when samples need **different** settings, or to add analyses/images beyond the profile defaults. Everything a run can express reduces to one canonical, list-based configuration that the three layers (**profile → CLI → JSON**) assemble:

```jsonc
{
  "platform": "10x_xenium", "out_dir": "...", "resources": { "n_jobs": 8, "threads": 16 },
  "samples": [ { "id": "s1", "in_dir": "...", "tissue": "lung" } ],  // same fields as sheet columns, plus "images"
  "exclude_feature_regex": "...",
  "include_feature_list": "...",                    // or "exclude_feature_list"; factor analyses only
  "ingest_exclude_feature_regex": "...",            // ingest_{include,exclude}_feature_{regex,list}: drops from the data
  "ingest": { "jitter_xy": 0.8, "units_per_um": 1 }, // random +/- um offset on X/Y, and the input's coordinate units (see above)
  "ficture_defaults": { "decode_scale": 2 },
  "cell_defaults":    { "min_cell_count": 20 },
  "segment_10x":   { "widths": "12,24" },           // or true: export the hexagon files as 10x MEX (see FICTURE mode)
  "tile_size": 5000, "tile_buffer": 1000,           // punkst tiling knobs; omit to keep run_ficture2_multi's 500 / 1000
  "ficture":       [ /* analyses: each is a de-novo train OR a projection */ ],
  "cell_analyses": [ /* see below */ ],
  "cell_lists":    { "clusters": "..." },           // ready-made --list-* file(s) for every cell analysis
  "images":        [ /* see Image Modalities */ ],
  "cartload":  { "use_pmpoint": true, "bin_count": 500 }
}
```

Publishing (annotation + S3 upload) is **not** part of this config — it is driven entirely by CLI flags (see [Overview → Publishing](./run_together.md#publishing)).

### Top-level keys

Dict values are **deep-merged** into the profile; other values replace it, except the three lists merged by `id` (next section).

| Key | Meaning |
|-----|---------|
| `platform`, `out_dir`, `out_root`, `saw`, `colname_cell` | As `--platform`, `--out-dir`, `--out-root`, `--saw`, `--colname-cell` (the CLI wins). |
| `samples` | Sample entries (below), added to any `--samples` rows and the single-sample CLI sample. |
| `exclude_feature_regex`, `include_feature_list`, `exclude_feature_list` | Analysis feature filters (the CLI wins). |
| `ingest_include_feature_regex`, `ingest_exclude_feature_regex`, `ingest_include_feature_list`, `ingest_exclude_feature_list` | Ingest feature filters (the CLI wins). |
| `cell_lists` | Run-wide ready-made `--list-*` files, role → path (a `--list-*` flag wins per role). |
| `ingest` | Ingest settings (below). |
| `roles` | Per-role auto-detection and column settings (below). |
| `cartload` | `use_pmpoint`, `sge_scale`, `bin_count` (`--bin-count` wins). |
| `resources` | `n_jobs`, `threads` — passed to every command; **wins over** `-j` / `--threads` (`make` itself still runs with the CLI `-j`). |
| `ficture_defaults` / `cell_defaults` | Defaults merged into every `ficture` / `cell_analyses` entry (per-entry keys win). |
| `ficture`, `cell_analyses`, `images` | Lists merged by `id` (below). |
| `segment_10x`, `tile_size`, `tile_buffer` | As `--segment-10x` / `--segment-width-10x`, `--tile-size`, `--tile-buffer`. |
| `image_defaults` | Fallback image settings; see [Image Modalities](./run_together_images.md#fields). |
| `image_transform` | `{ "column": ..., "flag": ... }`: the platform transform column of image rows (MERSCOPE: `merfish_csv` → `--micron2pixel-csv`). |
| `hne` | Settings for a sample's `hne` image (Visium HD profile: `georef_plain`, `um_per_pixel_json`, `um_per_pixel_key`, `out_id`). |
| `squares` | Visium HD square-bin imports, a list of `{ "subdir", "bin_size", "suffix" }`. |
| `cell_import` | Visium HD segmented-cell import, `{ "detect_dir", "suffix" }`. |
| `publish` | Accepted but not read (publishing is CLI-only). |

**Sample entries** accept `id`, `in_dir`, `in_prefix`, `raw_transcript`, the role keys and their aliases (`transcript`/`tsv`, `xy`/`cell_xy`, `boundaries`/`cell_boundary`/`cell_boundaries`, `clusters`, `mex`/`mex_dir` or `mex_bcd`+`mex_ftr`+`mex_mtx`, `cellxgene`/`cell_by_gene`, `cell_tsv`, `gef`, `cellbin_gef`), `dapi`, `hne`, `tissue`, and `images` — a list of image entries in the [config JSON form](./run_together_images.md#config-json-images) that take precedence over the profile's entries with the same `id`.

**`ingest`** keys:

| Key | Meaning |
|-----|---------|
| `method` | `sge_convert` (default), `reformat_cosmx`, or `stereoseq`. |
| `sge_platform` | `sge_convert --platform` (default: the profile name). |
| `units_per_um`, `jitter_xy` | As `--units-per-um`, `--jitter-xy` (the CLI wins). |
| `autodetect` | List of `{ "file", "flag" }`: raw transcript candidates under `in_dir`; the first that exists is used. |
| `inputs` | Flag → path under `in_dir` (`sge_convert`), or flag → glob (`reformat_cosmx`; `null` drops that input). |
| `input_roles` | Flag → role whose path is passed (Illumina: `--in-mex` ← `mex`). |
| `in_dir_flag` | Flag that receives `in_dir` itself (Seq-Scope: `--in-mex`); skipped when the sample has a `raw_transcript`. |
| `raw_input_flag` | Flag that receives `raw_transcript` (default `--in-csv`). |
| `produces` | `reformat_cosmx`: output suffix per role (`transcript`, `xy`, `boundaries`). |
| `csv_colnames` | `{ "x", "y", "feature", "count" }` raw transcript columns (the `--colname-transcript-*` flags). |
| `csv_colname_cell` | As `--colname-transcript-cell`. |
| `csv_colnames_others` | Extra raw columns carried into the transcript TSV (Xenium: `cell_id`, `z_location`, …). |
| `assign_cell_id` | `true`, or `{ "expand_um", "colname_x", "colname_y" }`: run `spatula tsv-add-cell-id` for samples with `boundaries` (`colname_x`/`colname_y` are the transcript's columns, default `X`/`Y`; `expand_um` assigns an outside-all-cells transcript to the nearest cell within that distance). |
| `cellxgene_blank_prefix` | `convert_cellxgene --blank-prefix` (its default: `Blank`). |
| `extra_flags` | Raw flags appended to the ingest command (for Stereo-seq, the bin1 `sge_convert`). Lists are replaced, not merged. |
| `gef_suffix` | Stereo-seq: only shown in error messages (the detected suffix is `roles.gef.suffix`). |
| `feature_file` | Stereo-seq: the feature file the cell-bin conversion checks against (default `feature.clean.tsv.gz`). |

**`roles`** — one entry per role: `file` (auto-detect at this path under `in_dir`) or `suffix` (auto-detect at `in_prefix` + suffix), and `colname_cell` / `colname_x` / `colname_y` (`xy`: passed to the cell decode; `boundaries`: read by `tsv-add-cell-id`). A `boundaries` entry may also set `format` (`geojson`: derive cell centroids from the polygons), `cell_id_format`, `cell_id_prop`, `units_key`, and `units_json` (scale JSON under `in_dir`) — see the [Visium HD profile](./platforms/visium_hd.md).

### Validation of config keys

A misspelled key would otherwise be merged in and never read, so `--config` is checked before planning; any of the following aborts the run with the list of recognized keys:

- an unrecognized **top-level** key (e.g. a top-level `csv_colnames`, which belongs under `ingest`);
- an unrecognized key in a **sample** entry, or in the closed-schema blocks **`ingest`**, **`cartload`**, **`resources`**, **`image_defaults`**, and `segment_10x`;
- a sample **role the platform cannot consume** (e.g. `cell_tsv` on `10x_xenium`) — usually a wrong `--platform` or a leftover column.

`--allow-unknown-config-keys` downgrades all of these to "ignore and continue". Open maps (`roles`, `hne`, `image_transform`, `cell_import`, `squares`, `publish`, per-entry list fields), `--platform-json`, and sample-sheet columns are not checked.

### List assembly rule (append-by-default, keyed by `id`)

For `ficture`, `cell_analyses`, and `images`, JSON entries are **merged into** the profile/CLI-derived list:

- `ficture`, `cell_analyses`: an entry with a **new `id`** is appended; one reusing an **existing `id`** is deep-merged into it.
- `images`: an entry **replaces every** profile entry with the same `id` (profiles may list several same-`id` fallbacks, e.g. Xenium's `dapi`); new ids are appended.
- `{ "replace": [ ... ] }` discards the base list entirely; `{ "add": [ ... ] }` is the same as a plain list.

This is what lets you set the base on the CLI and add only the extras in JSON:

```bash
cartloader run_together --platform 10x_xenium --in-dir IN --out-dir OUT \
    --width 12 --n-factor 24 --config extra.json
```
```jsonc
// extra.json — de-novo base comes from the CLI; these are ADDED
{
  "ficture":       [ { "id": "ref", "mode": "project", "model": "/models/ref.tsv", "width": 12 } ],
  "cell_analyses": [ { "id": "spatch", "uses": ["xy","boundaries","clusters","mex"], "model_id": "ref" } ],
  "images":        [ { "id": "cd3", "source": "cd3.ome.tif", "kind": "single", "color": "FF0000" } ]
}
```

### `ficture` analyses

Each entry is either de-novo or a projection:

```jsonc
{ "id": "denovo", "mode": "train",   "width": "12", "n_factor": "24,48,96" }
{ "id": "ref",    "mode": "project", "model": "/models/ref.tsv", "width": 12 }
```

| Key | Meaning |
|-----|---------|
| `id` | Entry id. A de-novo entry produces models `t<width>_f<n_factor>`; a projection's model id is its `id`. |
| `mode` | `train` (default) or `project`. |
| `width`, `n_factor` | Comma lists; de-novo trains every combination. On a projection, `n_factor` only informs which model is the default for cell analyses (largest wins). |
| `model` | Projection: path to the trained `*.model.tsv`. |
| `min_ct_per_unit_hexagon`, `min_ct_per_unit_train`, `decode_scale` | Per-analysis `run_ficture2_multi` parameters. |
| `alias` | A factor-label file copied beside the factor as `<id>-alias.tsv` (`_` → `-` in the id) and recorded under the factor's `alias` key. The entry's `id` must be the catalog's factor id, so this is for projection entries; the file is copied as given. |

`ficture_defaults` (e.g. `decode_scale`, `min_ct_per_unit_hexagon`, `min_ct_per_unit_train`) apply to **every** analysis, including projections; per-entry keys win. `cell_defaults` does the same for `cell_analyses` entries (`min_cell_count`, `min_feature_count`).

### `cell_analyses`

Cell-level decode is **platform-default and automatic**: an analysis runs whenever a sample provides its required roles. A sample contributes when it has **every** role in `uses` **and** (if present) **at least one** role in `any_uses`; roles in `optional_uses` are added to the decode when available.

```jsonc
{ "id": "spatch", "uses": ["xy", "boundaries", "clusters", "mex"], "model_id": "ref" }
```

| Key | Meaning |
|-----|---------|
| `id` | Analysis id; names the factor in the catalog. |
| `uses` / `any_uses` / `optional_uses` | Required (all) / alternative (at least one) / optional roles. The pseudo-role `tx_cell_id` matches a sample whose transcript already carries a `cell_id` column. |
| `require_all` | Run only when **every** sample in the group contributes. |
| `generic_cell` | `generic` platform: run from `--colname-cell` or from a `mex` present in every sample (not both); see [Generic](./platforms/generic.md). |
| `multi_import` | On joint runs, import per sample with this command instead of decoding jointly (Xenium: `import_xenium_cell`). |
| `model_id` | FICTURE model to decode against (default: the largest-factor model). |
| `min_cell_count`, `min_feature_count` | Per-analysis `run_ficture2_multi_cells` thresholds (else `cell_defaults`, else `--cell-min-*`). |
| `lists` | Ready-made `--list-*` files for this analysis, role → path (see below). |
| `extra_flags` | Raw flags appended to this analysis's `run_ficture2_multi_cells` call (e.g. `["--zero-based-clust-id"]`, `["--tsne"]`). |
| `name` | Display name in the catalogs (below). |
| `alias` | Companion cluster-label file (below). |

`any_uses` lets an analysis accept **alternative** cell-count sources — e.g. MERSCOPE runs cell analysis from **either** boundaries **or** a `cellxgene` MEX, and a **mixed joint run** (some samples with boundaries, some with a MEX) is resolved per sample and packaged from one call. See the [MERSCOPE page](./platforms/merscope.md) for the source-precedence rules.

A **`multi_import`** analysis (e.g. Xenium's `xeniumranger` → `import_xenium_cell`) relies on **sample-specific** cluster labels: on a single-sample run it goes through `run_ficture2_multi_cells` as usual; on a **joint run** it is imported **per sample** in the images stage instead (sheet-provided `xy`/`boundaries`/`clusters` paths are forwarded as `--csv-*` overrides). See the [Xenium page](./platforms/xenium.md).

TSNE manifolds are **not** generated by default (only UMAP is; nothing downstream reads TSNE, and it is slow on large datasets). To add them for an analysis, pass `"extra_flags": ["--tsne"]`.

`name` sets the factor's human-readable **`name:`** in `catalog.yaml` and `multi-catalog.yaml` — the label shown for the layer. The factor **`id`** is unchanged (it still names every file and asset key), so this is purely cosmetic:

```jsonc
{ "id": "published", "name": "Published cell types (Banovich 2025)", "uses": ["boundaries", "xy"] }
```

Without it, a cell analysis displays its bare id (`published`), and the multi-catalog's cell factors carry no `name` at all. The name is written after packaging, into the shared multi-catalog and every contributing sample's catalog. Plain text only — quotes, `$`, backticks, backslashes and newlines are rejected up front, since the name travels through a generated make recipe. Unrelated to `alias`, which points a factor at a companion factor-label TSV.

`alias` names a companion **cluster-label file** for the analysis — one label per cluster id:

```jsonc
{ "id": "published", "name": "Published cell types", "alias": "/work/clust/published.alias.tsv",
  "lists": { "clusters": "/work/clust/published.list.tsv" } }
```
```
cluster	annotation
1	EN-IT
2	EN-ET
3	IN
```

Two columns, cluster id then label, with an optional header; tab-, comma- or whitespace-separated (a label may contain spaces unless the file is comma-separated). The ids use **the same numbering as the analysis's cluster files**: 1-based for supplied cluster labels, unless the analysis passes `--zero-based-clust-id` in `extra_flags`, in which case they are 0-based. A Leiden analysis clusters on demand, so its alias ids are 0-based, matching the packaged output. At planning time the file is validated (integer ids at or above the base, no duplicates, non-empty labels) and rewritten to the packaged layout (`index<TAB>alias`, 0-based) as `<out_dir>/tsv/alias.<analysis_id>.tsv`. After packaging it is deployed as `<analysis_id>-alias.tsv` and recorded under the factor's `alias` key in the multi-catalog and every contributing sample's catalog. A factor that carries `alias` is skipped by the AI annotation stage (`--anno` / `--anno-deep`), so curated labels are not overwritten.

### Supplying your own `--list-*` files

By default `run_together` **derives** each `--list-*` file that [`run_ficture2_multi_cells`](./run_ficture2_multi_cells.md) consumes, writing `<out_dir>/tsv/in_<role>.<analysis_id>.tsv` from the samples' resolved roles. To supply one yourself instead — most often **externally assigned cell clusters**, since without `--list-cluster` the cells stage computes Leiden clusters on demand — name it per role, either run-wide or per analysis:

| role | flag it feeds | CLI flag | line format |
|---|---|---|---|
| `clusters` | `--list-cluster` | `--list-cluster` | `SAMPLE_ID<TAB>CLUSTER_FILE` |
| `xy` | `--list-xy` | `--list-xy` | `SAMPLE_ID<TAB>XY_FILE` |
| `boundaries` | `--list-boundaries` | `--list-boundaries` | `SAMPLE_ID<TAB>BOUNDARY_FILE` |
| `mex` | `--mex-list` | `--list-mex` | `SAMPLE_ID<TAB>MEX_DIR` (or a bcd/ftr/mtx triple) |
| `cell_tsv` | `--tsv-list` | `--list-cell-tsv` | `SAMPLE_ID<TAB>CELL_TSV` |

```bash
# run-wide, from the CLI (applies to every cell analysis)
cartloader run_together --platform 10x_xenium --samples samples.tsv --out-dir OUT \
    --list-cluster /work/clust/list.tsv
```
```jsonc
// run-wide, in JSON: same effect as the CLI flags (a CLI flag wins)
{ "cell_lists": { "clusters": "/work/clust/list.tsv" } }

// per analysis: overrides the run-wide default role by role
{ "cell_analyses": [ { "id": "cartloader",
                       "lists": { "clusters": "/work/clust/list.tsv" },
                       "extra_flags": ["--zero-based-clust-id"] } ] }
```

A named list is passed **verbatim** (no file is generated for that role) and:

- **satisfies that role's gating** — no sample has to carry the role on disk, and the role need not appear in the analysis's `uses` at all, so a cluster list can be attached to an analysis that would otherwise cluster on demand;
- is **validated at plan time** — unknown role name, missing file, empty file, and a first column matching none of the run's sample ids are hard errors; partial coverage and ids outside the run are warnings. The first column must use the run's sample ids (e.g. `rep1` for a lone `--out-dir` sample without `--id`);
- does **not** change `multi_import` routing — that path is per-sample by nature, so on a joint Xenium run attach the list to the jointly decoded `cartloader` analysis, and `xeniumranger` keeps its per-sample import.

Cluster ids are read as **1-based** and decremented; for 0-based labels add `"extra_flags": ["--zero-based-clust-id"]`.

---
## See also

- [`run_together` Overview](./run_together.md) — the pipeline, stages/resume, publishing, output.
- [Image Modalities](./run_together_images.md) — how to attach DAPI/H&E/protein images.
- [Supported Platforms](./supported_platforms.md) and the per-platform pages.
- Guides: [single sample](../vignettes/pipelines/run_together.md), [multiple samples](../vignettes/pipelines/run_together_multi.md).
