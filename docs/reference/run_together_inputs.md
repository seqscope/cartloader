# Specifying Inputs to `run_together`

`run_together` accepts inputs in **three ways**, from most convenient to most expressive. They compose — you can set the base run on the CLI and add only the extras in a config JSON.

| Mode | Flag | Use when |
|------|------|----------|
| **1. Single sample (direct CLI)** | `--in-dir` **or** explicit `--in-*` file flags | One sample; a standard platform directory, or a few individually-named files. |
| **2. Multi-sample (sample sheet)** | `--samples samples.tsv` | Several samples (a joint model by default). A wide TSV, one row per sample. |
| **3. Full config (JSON)** | `--config run.json` | Per-sample settings differ, or you need analyses/images beyond the profile defaults. |

Whichever you use, the **platform profile** (`--platform <name>`) supplies the defaults: which files to auto-detect, the default FICTURE analyses, and the image conventions. See [Supported Platforms](./supported_platforms.md) for what each profile expects and [per-platform pages](./platforms/xenium.md) for worked examples of all three modes.

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
    --in-transcript    /data/molecules.csv \
    --in-cell-xy       /data/cell_metadata.csv \
    --in-cell-boundary /data/cell_boundaries.csv \
    --id MYSAMPLE --out-dir OUT
```

| Flag | Sample-sheet column | Meaning |
|------|--------------------|---------|
| `--in-dir` | `in_dir` | Raw platform directory → every role auto-detected inside it. A platform may instead consume the directory whole (Seq-Scope: it *is* the MEX directory) |
| `--in-prefix` | `in_prefix` | Raw platform **path prefix** → roles auto-detected by suffix (Stereo-seq) |
| `--in-transcript` | `raw_transcript` | Raw transcript CSV/TSV to ingest (through `sge_convert`) |
| `--in-cell-xy` | `xy` | Cell centroids / metadata file |
| `--in-cell-boundary` | `boundaries` | Cell boundary polygons |
| `--in-cellxgene` | `cellxgene` | Cell×gene matrix CSV → converted to a MEX (drives cell clustering) |
| `--id` | `id` | Sample id. With `--out-dir`, defaults to **`rep1`** (the packaged dir/catalog id becomes `<out-dir basename>-rep1`); with `--out-root`, defaults to the `in_dir`/`in_prefix` basename |

`--in-dir` is just a one-row sample sheet with only `in_dir`; the explicit `--in-*` flags are a one-row sheet with those columns.

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

For **independent** per-sample models (each its own model), use `--out-root` instead of `--out-dir`; each row becomes its own `<out_root>/<id>/`.

### The sample sheet is a wide table of input roles

Columns map to per-sample **input roles**. Every column is optional except that each sample needs a transcript source (`in_dir`, `raw_transcript`, or `transcript`). Column names are case-sensitive; the listed **aliases** are accepted interchangeably.

| Column (aliases) | Role | Meaning |
|---|---|---|
| `id` | — | Sample identifier. Optional: defaults to the `in_dir`/`in_prefix` basename, or to `rep1` for a lone `--out-dir` sample with no named input. |
| `in_dir` | — | Raw platform directory → the profile **auto-detects** every role inside it. |
| `in_prefix` | — | Raw platform path prefix → the profile auto-detects roles by **suffix** (`{prefix}.tissue.gef`, …). For platforms whose files share a name rather than a directory; see [BGI Stereo-seq](./platforms/stereoseq.md). |
| `raw_transcript` | — | Explicit path to a **raw** transcript file to ingest (e.g. MERSCOPE `detected_transcripts.csv`). Runs through `sge_convert`. Sheet equivalent of `--in-transcript`. |
| `transcript` (`tsv`) | transcript | A **pre-converted** `transcripts.tsv.gz` (the *output* of `sge_convert`) → **skips ingest**. Not for raw CSVs. |
| `xy` (`cell_xy`) | xy | Cell centroids file. |
| `boundaries` (`cell_boundary`, `cell_boundaries`) | boundaries | Cell boundary polygons. |
| `clusters` | clusters | External cluster labels. |
| `mex` (`mex_dir`) | mex | MEX directory. Or give the explicit triple `mex_bcd` / `mex_ftr` / `mex_mtx` when one directory does not apply. |
| `cellxgene` (`cell_by_gene`) | cellxgene | Cell×gene matrix CSV (e.g. MERSCOPE `cell_by_gene.csv`) → converted to a MEX that drives cell clustering (works without boundaries). |
| `gef` / `cellbin_gef` | gef, cellbin_gef | Stereo-seq binary GEFs, when they do not match `in_prefix` + the standard suffix. |
| `cell_tsv` | cell_tsv | A standalone pixel TSV (`X`, `Y`, gene, count, cell_id) that supplies cell counts on its own, for platforms whose cell assignment cannot be carried on the transcript (Stereo-seq cell bins). |
| `hne` | — | H&E image (Visium HD adds the layer automatically). |

Images are **not** sample-sheet columns (except `hne`) — they are a separate concern; see [Image Modalities](./run_together_images.md).

**Rules:** an explicit column **overrides** auto-detection for that role; the cell values `` (empty), `-`, `.`, and `NA` all mean *unset*. Role columns are resolved relative to `in_dir` when relative, but `raw_transcript` is not — give it an absolute path (or one relative to the working directory).

### The three transcript sources

Each sample gets its transcripts from exactly one of:

- **`in_dir`** — the standard path: point at the raw platform folder and the profile finds `detected_transcripts.csv[.gz]` (and every other role) inside it.
- **`raw_transcript`** — a raw CSV/TSV that still needs ingesting, when it is arbitrarily named or not laid out as a standard `in_dir`.
- **`transcript`** — an already-ingested TSV (`transcripts.tsv.gz`); ingest is skipped and it feeds FICTURE directly.

```
# in_dir (auto-detect), a pre-converted TSV, and a raw CSV named explicitly:
id    in_dir            transcript                    raw_transcript                 boundaries
s1    /data/s1
s2                      /data/s2/transcripts.tsv.gz
s3                                                    /data/s3/detected.csv.gz       /data/s3/bounds.csv.gz
```

### Column-name overrides

When input **columns** are non-standard, name them (CLI or the profile/`--config`):

- `--colname-transcript-x/-y/-feature/-count` — the raw transcript's coordinate/gene/count columns.
- `--colname-xy-cell/-x/-y` — the xy file's columns (use `--colname-xy-cell ''` for an unnamed index column).
- `--colname-boundary-cell/-x/-y` — the boundary file's columns.

!!! note "Transcript that already carries a `cell_id` column (rare)"
    If your raw transcript CSV already has a per-molecule cell-id column, name it with
    `--colname-transcript-cell <name>` (or `ingest.csv_colname_cell` in `--config`).
    The column is carried through ingest to transcript column 5 and used directly for
    cell analysis, skipping `spatula tsv-add-cell-id`. A `cellxgene` MEX, if also
    present, still drives the clustering.

---
## FICTURE mode (de-novo vs. projection)

Independent of how inputs are specified, Tier-1 selects the base FICTURE work. Exactly one mode is the base (a config JSON can add more analyses on top).

=== "De-novo training (default)"

    Train new LDA models. `--width` and `--n-factor` accept comma lists → the cross-product is trained (**multiple widths supported**).

    ```bash
    cartloader run_together --platform 10x_xenium --in-dir IN --out-dir OUT \
        --width 12,18 --n-factor 24,48
    ```

=== "Projection-only"

    Reuse **already-trained** models — no LDA training runs. Point `--project-models` at one or more **existing FICTURE directories**; `run_together` reads each `ficture.params.json` and re-projects every model it lists onto the current data.

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

    Applies to **every** platform, not just Seq-Scope. Cell analyses are skipped too (they decode against a model). Cannot be combined with `--n-factor` or `--project-models`. The hexagon files are still built alongside the tiles, so adding factors later re-uses them instead of re-tiling — rerun without `--no-ficture` in the same `--out-dir`.

**Common decode overrides** (else profile / built-in default): `--exclude-feature-regex`, `--min-ct-per-unit-hexagon` (default `50`), and `--always-single-molecule` / `--never-single-molecule` (default: single-molecule **ON** for pixel FICTURE, **OFF** for cell decode). An explicit CLI flag wins over a `--config`/profile value, which wins over the built-in default.

---
## Mode 3 — Full config (JSON)

Escalate to `--config run.json` when samples need **different** settings, or to add analyses/images beyond the profile defaults. Everything a run can express reduces to one canonical, list-based configuration that the three layers (**profile → CLI → JSON**) assemble:

```jsonc
{
  "platform": "10x_xenium", "out_dir": "...", "resources": { "n_jobs": 8, "threads": 16 },
  "samples": [ { "id": "s1", "in_dir": "..." } ],   // same fields as sheet columns: in_dir, raw_transcript, transcript, xy, boundaries, clusters, mex, cellxgene
  "exclude_feature_regex": "...",
  "ficture_defaults": { "decode_scale": 2 },
  "ficture":       [ /* analyses: each is a de-novo train OR a projection */ ],
  "cell_analyses": [ /* {id, uses:[roles], any_uses?:[roles], optional_uses?:[roles], model_id?} */ ],
  "images":        [ /* see Image Modalities */ ],
  "cartload":  { "use_pmpoint": true, "bin_count": 500 }
}
```

Publishing (annotation + S3 upload) is **not** part of this config — it is driven entirely by CLI flags (see [Overview → Publishing](./run_together.md#publishing)).

### List assembly rule (append-by-default, keyed by `id`)

For `ficture`, `cell_analyses`, and `images`, JSON entries are **merged into** the profile/CLI-derived list:

- entry with a **new `id`** → appended;
- entry reusing an **existing `id`** → deep-merged (override);
- to discard the base list entirely, write the section as `{ "replace": [ ... ] }`.

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

`ficture_defaults` (per-analysis decode params like `decode_scale`) apply to **every** analysis, including projections; per-entry keys win.

### `cell_analyses`

Cell-level decode is **platform-default and automatic**: an analysis runs whenever a sample provides its required roles. A sample contributes when it has **every** role in `uses` **and** (if present) **at least one** role in `any_uses`; roles in `optional_uses` are added to the decode when available. `model_id` picks which FICTURE model decodes the cells (default: the largest-factor model).

```jsonc
{ "id": "spatch", "uses": ["xy", "boundaries", "clusters", "mex"], "model_id": "ref" }
```

`any_uses` lets an analysis accept **alternative** cell-count sources — e.g. MERSCOPE runs cell analysis from **either** boundaries **or** a `cellxgene` MEX, and a **mixed joint run** (some samples with boundaries, some with a MEX) is resolved per sample and packaged from one call. See the [MERSCOPE page](./platforms/merscope.md) for the source-precedence rules.

An analysis may carry a **`multi_import`** command (e.g. Xenium's `xeniumranger` → `import_xenium_cell`). Such an analysis relies on **sample-specific** cluster labels: on a single-sample run it goes through `run_ficture2_multi_cells` as usual; on a **joint run** it is imported **per sample** instead (sheet-provided `xy`/`boundaries`/`clusters` paths are forwarded as `--csv-*` overrides). See the [Xenium page](./platforms/xenium.md).

---
## See also

- [`run_together` Overview](./run_together.md) — the pipeline, stages/resume, publishing, output.
- [Image Modalities](./run_together_images.md) — how to attach DAPI/H&E/protein images.
- [Supported Platforms](./supported_platforms.md) and the per-platform pages.
- Tutorials: [single-sample](../vignettes/pipelines/run_together.md), [multi-sample](../vignettes/pipelines/run_together_multi.md).
