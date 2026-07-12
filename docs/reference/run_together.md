# Unified Pipeline Orchestrator (`run_together`)

## Overview

`run_together` runs a complete CartoScope pipeline — **ingest → FICTURE → cell decode → asset packaging → image import → (optional) publish** — across many spatial platforms, for one or many samples, from a single command. It assembles the whole run as **one Makefile** (`run_together.mk`) and executes it, so the pipeline is **resumable**, **parallel**, and trains **one joint FICTURE model** across a multi-sample run.

The design has one guiding principle:

> **Convenience for the default case; full control in JSON.**

A standard run needs only a few flags. Everything a run *can* express lives in a single canonical, list-based configuration that three layers assemble:

1. **Profile** — per-platform defaults (auto-detected inputs, default analyses, image conventions).
2. **Tier-1 CLI** — the common knobs (`--in-dir`/`--samples`, `--width`/`--n-factor` or `--project-models`, `--out-dir`).
3. **Tier-2 JSON** (`--config`) — augments or fully specifies anything.

The layers **compose**: you can set the base run on the CLI and add only the extra analyses/images in JSON.

---
## Requirements

- `make` on the `PATH`.
- Tools used by the stages: `spatula`, `punkst` (FICTURE2), `pigz`, `sort`, `python3`, `go-pmtiles`, `gdal`, `tippecanoe`, `parquet-tools`, `jq` (Visium HD H&E), and `aws` (publish). `run_together` delegates to the CartLoader modules, which find `spatula`/`punkst` as built repo submodules or on `PATH`.

---
## Input: single vs. multi-sample

One symmetric choice; everything else is identical:

| | Flag |
|---|---|
| Single sample | `--in-dir DIR` |
| Multi-sample (joint model) | `--samples samples.tsv` |

Multi-sample defaults to a **joint model** (all samples share one `out_dir`). Independent per-sample models are an explicit opt-in with `--out-root` (each sample gets its own directory and model).

### The sample sheet is a wide table of input roles

Columns map to per-sample **input roles**. All are optional except `id`:

| Column | Meaning |
|--------|---------|
| `id` | Sample identifier (required) |
| `in_dir` | Raw platform directory → the profile **auto-detects** every role inside it |
| `transcript` | A pre-converted `transcripts.tsv.gz` → **skips ingest** for that sample |
| `xy` | Cell centroids file |
| `boundaries` | Cell boundaries file |
| `clusters` | External cluster labels |
| `mex` | MEX directory |

```
id    in_dir            transcript                    boundaries
s1    /data/s1
s2                      /data/s2/transcripts.tsv.gz   /data/s2/bounds.csv.gz
```

Rules: an explicit column **overrides** auto-detection for that role; if `transcript` is given, `sge_convert` is skipped and the TSV feeds FICTURE directly; `--in-dir` is just a one-row sheet with only `in_dir`.

---
## The two FICTURE modes

Tier-1 selects the base FICTURE work. Exactly one mode is the base (JSON can add more analyses on top).

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

---
## The canonical configuration

Everything reduces to these sections. The three list sections are **assembled across the profile, CLI, and JSON layers**.

```jsonc
{
  "platform": "10x_xenium", "out_dir": "...", "resources": { "n_jobs": 8, "threads": 16 },
  "samples": [ { "id": "s1", "in_dir": "..." } ],
  "exclude_feature_regex": "...",
  "ficture_defaults": { "decode_scale": 2 },
  "ficture":       [ /* analyses: each is a de-novo train OR a projection */ ],
  "cell_analyses": [ /* {id, uses:[roles], model_id?} */ ],
  "images":        [ /* {id, source|match, kind, color, convert} */ ],
  "cartload":  { "use_pmpoint": true, "bin_count": 500 }
}
```

Publishing (annotation + S3 upload) is **not** part of this config — it is driven entirely by CLI flags (see [Publish](#publishing)).

### List assembly rule (append-by-default, keyed by `id`)

For `ficture`, `cell_analyses`, and `images`, JSON entries are **merged into** the profile/CLI-derived list:

- entry with a **new `id`** → appended;
- entry reusing an **existing `id`** → deep-merged (override);
- to discard the base list entirely, write the section as `{ "replace": [ ... ] }`.

This is what lets you set the base on the CLI and add only the extras in JSON.

### `ficture` analyses

Each entry is either de-novo or a projection:

```jsonc
{ "id": "denovo", "mode": "train",   "width": "12", "n_factor": "24,48,96" }
{ "id": "ref",    "mode": "project", "model": "/models/ref.tsv", "width": 12 }
```

`ficture_defaults` (per-analysis decode params like `decode_scale`) apply to **every** analysis, including projections; per-entry keys win. The common `min_ct_per_unit_hexagon` and single-molecule behavior are instead controlled globally by the CLI flags below (`--min-ct-per-unit-hexagon`, `--always/never-single-molecule`).

### `cell_analyses`

Cell-level decode is **platform-default and automatic**: an analysis runs whenever every role in its `uses` list is available for a sample. Built-in profiles wire the standard ones (Xenium → `cartloader` + `xeniumranger`). Custom analyses (e.g. externally-clustered cells) are added in JSON:

```jsonc
{ "id": "spatch", "uses": ["xy", "boundaries", "clusters", "mex"], "model_id": "ref" }
```

`model_id` picks which FICTURE model decodes the cells (default: the largest-factor model). The role paths come from each sample's resolved roles (sheet columns / auto-detection).

### `images`

Profiles auto-detect standard modalities; per-sample/JSON entries append. Each declares a **kind**:

| `kind` | Behavior | Typical use |
|--------|----------|-------------|
| `single` | grayscale → colorize with `color` | DAPI, protein/RNA stains (OME-TIFF) |
| `rgb` | passthrough (no colorize), georeference | H&E histology |
| `prebuilt` | copy an existing `.pmtiles` into the catalog | pre-rendered layers |

```jsonc
{ "id": "dapi", "match": "morphology_focus/..._0000.ome.tif", "kind": "single", "color": "0F73E6", "convert": "ome2png" }
{ "id": "hne",  "source": "he.tif", "kind": "rgb", "convert": "png2pmtiles" }
{ "id": "boundary", "source": "prebuilt/boundary.pmtiles", "kind": "prebuilt" }
```

`source` is an explicit path (or `in_dir`-relative); `match` is an `in_dir`-relative pattern used for auto-detection. Same-`id` fallback entries are allowed — the first whose file exists is used.

---
## Mixing CLI + JSON (the common pattern)

Set the base run on the CLI; add only what's extra in JSON:

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

---
## Stage selection & resume

- **Resume:** re-run the same command (or `make -f OUT/run_together.mk -j N`); completed stages are skipped via flag files.
- `--only ingest,ficture` / `--skip images` — run a subset (stages: `ingest,ficture,cells,cartload,images,anno,upload`); excluded upstream stages are assumed done (prereqs are pruned so `make` won't error).
- `--restart` — rebuild everything (`make -B`).
- `--dry-run` — write the Makefile and print commands (`make -n`) without executing.

**Packaging bundles everything produced** — every FICTURE pixel decode plus every cell analysis that ran.

---
## Publishing

Publishing is **opt-in** and **entirely CLI-driven** (no config block). It has two independent actions — enable either or both:

- **`--anno`** — AI-annotate each packaged sample directory. **Requires `--tissue` and `--organism`** (no defaults); `--anno-api-type` (default `umgpt`), `--anno-model` (default `claude-opus-4-7`), and `--anno-threads` (default `10`) are overridable.
- **`--s3-upload`** — upload each self-contained sample directory to S3.

When both run, `upload` waits on `anno` (which edits `catalog.yaml`).

**S3 destination:** `<s3-prefix>/batch=<batch>/<collection>/<dir-id>/`, where `<dir-id>` is the sample's output directory name (`<sample_id>` for a single run, `<multi_id>-<sample_id>` for a joint run).

- `--s3-prefix` — default `s3://cartostore/data`
- `--batch` — default: the **current** `YYYY_MM`
- `--collection` — default: the **out-dir basename** (the run id)
- `--aws-profile` — AWS CLI profile (default `cartostore`); `--aws` — path to the `aws` binary (default `aws`)

```bash
cartloader run_together --platform 10x_xenium --samples samples.tsv --out-dir OUT \
    --width 12 --n-factor 24 \
    --anno --tissue "Kidney" --organism human \
    --s3-upload --collection my-collection
```

---
## How samples are packaged

- **Single sample** → one `run_cartload2` call, output at `cartl/<id>/`.
- **Joint multi-sample run** (several samples sharing one `--out-dir`) → a single [`run_cartload2_multi`](./run_cartload2_multi.md) call that packages every sample in parallel. Each sample is written to a **self-contained** `cartl/<multi_id>-<sample_id>/` directory (with `<multi_id>` defaulting to the `--out-dir` basename), and a **`cartl/multi-catalog.yaml`** is written that links every per-sample `catalog.yaml` and copies the shared factor files (`post`/`rgb`/`de`/`info`/`umap`) into the `cartl/` root as a unified `factors:` map, so `cartl/` uploads to S3 as one deployable unit.

This mirrors the FICTURE manifests: `run_ficture2_multi` writes a shared [`ficture.multi.params.json`](./run_ficture2_multi.md#multi-sample-manifest) that `run_cartload2_multi` reads to discover samples and shared assets.

---
## Command-line parameters

**Run:** `--dry-run`, `--restart`, `-j/--n-jobs`, `--threads`, `--makefn`, `--only`, `--skip`.

**Input/output:** `--platform`, `--in-dir`, `--samples`, `--out-dir`, `--out-root`, `--id`, `--config`, `--platform-json` (external profile override).

**Publish:** `--anno`, `--s3-upload`, `--tissue`, `--organism`, `--anno-api-type`, `--anno-model`, `--anno-threads`, `--collection`, `--batch`, `--s3-prefix`, `--aws-profile`, `--aws` (see [Publish](#publishing)).

**FICTURE mode:** `--width`, `--n-factor` (de-novo); `--project-models` (projection-only — existing FICTURE dir(s), comma-separated).

**Common decode overrides** (else profile / built-in default):

- `--exclude-feature-regex` — regex of features to exclude. Default: the profile's, else `^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated|System|Gm[0-9]|MT-|mt-|Rps|Rpl|NCS-|NCP-)`.
- `--min-ct-per-unit-hexagon` — minimum count per hexagon for FICTURE. Default: `50`.
- `--always-single-molecule` / `--never-single-molecule` — force single-molecule ON (or OFF) for **both** pixel FICTURE and cell decode. Default (neither flag): **ON for pixel FICTURE, OFF for cell decode**.

Precedence for the regex and min-count: an explicit CLI flag wins over a `--config`/profile value, which wins over the built-in default.

---
## Output

```
out_dir/
├── run_together.mk               # the generated pipeline
├── run_together.resolved.json    # fully-assembled config (provenance)
├── mk/                           # per-stage flag files
├── tsv/                          # per-sample SGE + role list files (in_list.tsv, in_<role>.<analysis>.tsv)
├── fic/                          # FICTURE results; per sample under fic/samples/<id>/
│   └── ficture.multi.params.json     # shared multi-sample manifest (joint runs)
└── cartl/                        # ← deploy this directory
    ├── multi-catalog.yaml            # joint runs only: links per-sample catalogs + shared assets
    ├── <multi_id>-<sample_id>/       # joint run: one self-contained dir per sample
    │   └── catalog.yaml + PMTiles
    └── <id>/                         # single-sample run: catalog.yaml + PMTiles
```

For a joint run, upload the whole `cartl/` directory: `multi-catalog.yaml` plus each self-contained `<multi_id>-<sample_id>/`. For `--out-root` batches (independent per-sample models), this layout is created under `<out_root>/<id>/`, and the master Makefile is written at `<out_root>/`.

See also: [Supported Platforms & Inputs](./supported_platforms.md), and the [single-sample](../vignettes/pipelines/run_together.md) and [multi-sample](../vignettes/pipelines/run_together_multi.md) tutorials.
