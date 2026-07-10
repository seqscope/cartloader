# Unified Pipeline Orchestrator (`run_together`)

## Overview

`run_together` runs a complete CartoScope pipeline — **ingest → FICTURE → cell decode → asset packaging → image import → (optional) publish** — across many spatial platforms (10x Xenium, 10x Visium HD, and more), for one or many samples, from a single command.

Instead of executing each step directly, `run_together` **emits one master `Makefile`** (`run_together.mk`) whose targets are per-stage flag files wired together with their true dependencies, then runs `make`. This gives you three things for free:

- **Robustness / resume** — a failed or interrupted run continues from where it stopped; completed stages are not repeated.
- **Parallelism** — independent samples and stages run concurrently under `make -j`.
- **A joint FICTURE model** — for multi-sample runs, one model is trained on all samples together, then reused per sample.

The boilerplate that differs between platforms (which input files to read, exclusion regexes, coordinate scaling, morphology image naming, packaging flags) is captured in a **platform profile**, so simple runs need only a few arguments.

For a two-sample run, the dependency graph fans in on a joint model, then fans back out per sample:

```text
ingest (sample 1) ─┐
                   ├─► ficture (joint model) ─┬─► cartload 1 ─► images 1 ─► publish 1
ingest (sample 2) ─┘                          └─► cartload 2 ─► images 2 ─► publish 2
```

---
## Requirements

- Input data from a supported platform (see [Platform Profiles](#platform-profiles)).
- `make` on the `PATH`.
- The tools used by the underlying steps, depending on which stages run: `spatula`, `punkst`, `pigz`/`gzip`, `sort`, `python3`, `go-pmtiles`, `gdal`, `tippecanoe`, `parquet-tools`, `jq` (Visium HD H&E), and `aws` (publish).

---
## Pipeline stages

`run_together` builds up to six stages. Each stage maps to an existing CartLoader module.

| Stage | Module | Purpose |
|-------|--------|---------|
| `ingest` | [`sge_convert`](./sge_convert.md) | Convert raw platform output to a unified transcript TSV |
| `ficture` | [`run_ficture2_multi`](./run_ficture2_multi.md) | Train a joint FICTURE model (or project a pretrained model) |
| `cells` | `run_ficture2_multi_cells` | Segmentation/cluster-based decode (only if the profile defines cells) |
| `cartload` | [`run_cartload2`](./run_cartload2.md) | Package SGE + FICTURE results into PMTiles and `catalog.yaml` |
| `images` | [`import_image`](./import_image.md), [`import_square`](./import_square.md), [`import_cell`](./import_cell.md), `image_png2pmtiles` | Import morphology / H&E / square / cell layers and append them to the catalog |
| `publish` | annotation + `aws s3 cp` | **Opt-in.** AI annotation and upload of the catalog to S3 |

The stage names above are exactly what you pass to `--only` / `--skip`.

---
## Platform Profiles

A profile is a JSON file that tells `run_together` how to read a platform's output and how to package it. Built-in profiles live in `assets/run_together_profiles/`.

!!! info "Available built-in profiles"
    | `--platform` | Expected input directory | Notes |
    |--------------|--------------------------|-------|
    | `10x_xenium` | Xenium Ranger output | morphology images + cell/cluster decode wired automatically |
    | `10x_visium_hd` | Space Ranger output | 2 µm binned SGE, square/cell imports, optional H&E |

    Additional platforms (MERSCOPE, CosMx, Stereo-seq, SeqScope) are supported by `sge_convert` today and profiles for them are planned; until then use a custom `--profile` file (see [Overriding a profile](#overriding-a-profile)).

### What each profile expects on disk

The convenience of the simple invocation comes from these fixed-path assumptions. When testing on real data, verify your input directory matches.

=== "`10x_xenium`"

    Under `--in-dir` (the Xenium Ranger output directory):

    | Purpose | Path (first match wins) |
    |---------|-------------------------|
    | Transcripts | `transcripts.csv.gz`, or `transcripts.parquet`, or `transcripts/transcripts.parquet` |
    | Cell boundaries | `cell_boundaries.csv.gz` |
    | Cell centroids | `cells.csv.gz` (columns `x_centroid`, `y_centroid`) |
    | Cluster labels | `analysis/clustering/gene_expression_graphclust/clusters.csv` (imported under the `xeniumranger` prefix) |
    | Morphology images | `morphology_focus/morphology_focus_000{0,1,2,3}.ome.tif` → `dapi`/`boundary`/`rna`/`protein`; or single `morphology_focus.ome.tif` / `morphology.ome.tif` → `dapi` |

    Missing optional files are simply skipped.

=== "`10x_visium_hd`"

    Under `--in-dir` (the Space Ranger `outs/` directory):

    | Purpose | Path |
    |---------|------|
    | Scale factors | `binned_outputs/square_002um/spatial/scalefactors_json.json` |
    | Count matrix | `binned_outputs/square_002um/filtered_feature_bc_matrix/` |
    | Bin positions | `binned_outputs/square_002um/spatial/tissue_positions.parquet` |
    | Square layers | `binned_outputs/square_008um`, `binned_outputs/square_016um` (imported if present) |
    | Segmented cells | `segmented_outputs/` (imported if present) |
    | H&E image | Provided per sample via the `hne` field (µm/pixel read from the 2 µm `scalefactors_json.json`) |

    Coordinates are scaled by 2 (`--scale-xy 2.0` at ingest, `--sge-scale 2` at packaging).

---
## Usage

`run_together` has three input tiers. Pick the simplest that fits.

!!! warning "Replace placeholders"
    Replace example paths and IDs before running. All examples use `--dry-run`; **remove it to execute**.

### Tier 1 — single sample from the command line

The minimal case. Everything not given comes from the platform profile.

=== "Xenium"

    ```bash
    cartloader run_together \
      --platform 10x_xenium \
      --in-dir  /path/to/xenium_ranger/outs \
      --out-dir /path/to/out/batch/collection/my-xenium-id \
      --n-factor 12,24,48 \
      -j 4 --threads 8 \
      --dry-run
    ```

=== "Visium HD"

    ```bash
    cartloader run_together \
      --platform 10x_visium_hd \
      --in-dir  /path/to/spaceranger/outs \
      --out-dir /path/to/out/batch/collection/my-visiumhd-id \
      --n-factor 24,48,96 \
      -j 4 --threads 8 \
      --dry-run
    ```

=== "Project a pretrained model"

    Supplying `--pretrained-model` switches FICTURE from de-novo training to projection (`--n-factor` is ignored).

    ```bash
    cartloader run_together \
      --platform 10x_xenium \
      --in-dir  /path/to/xenium_ranger/outs \
      --out-dir /path/to/out/batch/collection/my-xenium-id \
      --pretrained-model /path/to/model.tsv \
      -j 4 --dry-run
    ```

!!! tip "ID / collection / batch inference"
    For a single-sample run with no `--id`, the sample ID is the basename of `--out-dir`. If `--out-dir` follows the `.../<batch>/<collection>/<id>` convention, those parts are recognized as such.

### Tier 2 — JSON configuration (multi-sample, joint model)

Use a JSON config when you have several samples that should share one FICTURE model, or want to set defaults once and override per sample. See the [Configuration reference](#configuration-reference) for all keys.

```json
{
  "platform": "10x_visium_hd",
  "out_dir": "/path/to/out/2026_07/vhd-prostate/poc",
  "resources": { "n_jobs": 10, "threads": 24 },
  "defaults": {
    "ficture":  { "n_factor": "24,48,96" }
  },
  "samples": [
    { "id": "38088", "in_dir": "/data/SI_38088/outs", "hne": "/data/HE_38088.tif" },
    { "id": "39685", "in_dir": "/data/SI_39685/outs", "hne": "/data/HE_39685.tif" }
  ]
}
```

```bash
cartloader run_together --config run.json -j 10 --dry-run
```

All samples in a single `out_dir` share one joint FICTURE model; each sample is then packaged independently.

### Tier 3 — TSV sample sheet (batches)

For many independent samples, list them in a tab-separated sheet and give an output root. Each row becomes its own output directory under `--out-root` with its own (single-sample) model.

**`samples.tsv`**

```tsv
id      platform    in_dir                     pretrained_model    hne
repA    10x_xenium  /data/xen/A                -                   -
repB    10x_xenium  /data/xen/B                -                   -
sc1     10x_xenium  /data/xen/C                /models/ref.tsv     -
```

```bash
cartloader run_together \
  --sheet samples.tsv \
  --config defaults.json \
  --out-root /path/to/out/2026_07 \
  -j 6 --dry-run
```

- A cell value of `-` (or blank) means "use the default / profile value".
- Recognized columns: `id`, `platform`, `in_dir`, `hne`, `pretrained_model`, `model_id`, `n_factor`, `width` (the last four are folded into each sample's FICTURE settings). The flat scalar columns are a convenience; richer per-sample overrides belong in a `--config` file.
- `--config` is optional and supplies shared defaults for every row.

---
## Overriding a profile

Profile values can be overridden two ways; both are deep-merged, so you only specify the keys you change.

=== "External profile file"

    Point `--profile` at a JSON file that overrides the built-in profile for the platform.

    ```bash
    cartloader run_together --platform 10x_xenium \
      --profile my_xenium_overrides.json \
      --in-dir ... --out-dir ... --n-factor 12,24 --dry-run
    ```

=== "Inline in the config"

    A top-level `profile` block in the JSON config is merged over the built-in profile for all samples.

    ```json
    {
      "platform": "10x_xenium",
      "out_dir": "...",
      "profile": {
        "ficture": { "width": "18" },
        "exclude_feature_regex": "^(Custom|Neg)"
      },
      "samples": [ { "id": "id1", "in_dir": "..." } ]
    }
    ```

**Merge order (later wins):** built-in profile → `--profile` file → config `profile` block → config `defaults` → per-sample overrides.

---
## Selecting stages (resume / partial runs)

Because everything is a `make` target keyed on flag files, you can re-run subsets safely.

- `--only ingest,ficture` — run only these stages.
- `--skip images,publish` — run everything except these.
- `--restart` — ignore existing outputs and rebuild all selected stages (`make -B`).

When an upstream stage is excluded, its dependency is dropped from downstream targets — `run_together` assumes that stage's outputs already exist on disk, so `make` will not error looking for a rule to build them. This is what makes "resume from cartload" work:

```bash
# FICTURE already finished; just (re)build packaging and images
cartloader run_together --config run.json --skip ingest,ficture,cells
```

---
## Publishing (opt-in)

The `publish` stage (AI annotation + S3 upload) runs **only** when both a `publish` block is present in the config **and** the `--publish` flag is passed. It is never triggered by a Tier-1 command.

```json
"publish": {
  "collection": "vhd-prostate",
  "batch": "2026_07",
  "annotate": { "tissue": "Prostate cancer", "organism": "human",
                "api_type": "umgpt", "model": "claude-opus-4-7", "threads": 10 },
  "upload":   { "s3_prefix": "s3://cartostore/data", "profile": "cartostore" }
}
```

```bash
cartloader run_together --config run.json --publish -j 10
```

Assets are uploaded to `<s3_prefix>/batch=<batch>/<collection>/<id>/`.

---
## Configuration reference

Top-level keys of the JSON config:

| Key | Type | Description |
|-----|------|-------------|
| `platform` | string | Default platform for all samples (overridable per sample). |
| `out_dir` | string | Single output directory; all samples share one joint model. |
| `out_root` | string | Output root for batches; each sample gets `<out_root>/<id>`. Set one of `out_dir`/`out_root`. |
| `resources` | object | `{ "n_jobs": int, "threads": int }`. Defaults come from CLI `-j` / `--threads`. |
| `profile` | object | Inline profile override, deep-merged over the built-in profile. |
| `defaults` | object | Default `ficture` / `cartload` / `annotate` sub-objects applied to every sample. |
| `samples` | array | One object per sample (see below). Required. |
| `publish` | object | Publish settings (see [Publishing](#publishing-opt-in)). |

**Sample object**

| Key | Description |
|-----|-------------|
| `id` | Sample identifier (defaults to `basename(in_dir)`; or `basename(out_dir)` for a single-sample `out_dir` run). |
| `in_dir` | Input directory for the sample. **Required.** |
| `platform` | Overrides the top-level platform for this sample (enables mixed-platform batches). |
| `hne` | Path to an H&E image (Visium HD). |
| `ficture`, `cartload`, `annotate` | Per-sample override sub-objects. |
| `exclude_feature_regex` | Per-sample feature-exclusion regex. |
| `pretrained_model`, `model_id`, `n_factor`, `width` | Convenience scalars folded into this sample's `ficture` settings. |

**`ficture` sub-object**: `width`, `n_factor`, `min_ct_per_unit_hexagon`, `decode_scale`, `single_molecule`, `pretrained_model`, `model_id`.

**`cartload` sub-object**: `use_pmpoint`, `bin_count`, `sge_scale`.

---
## Command-line parameters

### Run options

- `--dry-run`: Write the Makefile and print the commands (`make -n`) without executing.
- `--restart`: Ignore existing outputs and rebuild (`make -B`).
- `-j`, `--n-jobs`: Parallel jobs for `make` and passed to sub-commands (default: 1).
- `--threads`: Threads per job (default: 4).
- `--makefn`: Master Makefile name (default: `run_together.mk`).
- `--only`: Comma-separated stages to run exclusively.
- `--skip`: Comma-separated stages to skip.
- `--publish`: Enable the opt-in publish stage (requires a `publish` config block).

### Input / output

- `--config`: JSON run configuration (Tier 2).
- `--sheet`: TSV sample sheet (Tier 3).
- `--profile`: External JSON profile overriding the built-in one.
- `--platform`: Platform preset (Tier 1), e.g. `10x_xenium`, `10x_visium_hd`.
- `--in-dir`: Input directory for a single-sample Tier-1 run.
- `--out-dir`: Output directory (single run / shared joint model).
- `--out-root`: Output root under which each `--sheet` sample gets its own directory.
- `--id`: Sample ID for a single-sample Tier-1 run.

### Common FICTURE overrides

- `--n-factor`: Comma-separated factor counts (overrides the profile).
- `--width`: Hexagon width in µm (overrides the profile).
- `--pretrained-model`: Pretrained model TSV → projection instead of de-novo training.

---
## Output

Under each `out_dir`:

```
out_dir/
├── run_together.mk               # master Makefile (the pipeline)
├── run_together.resolved.json    # fully-resolved per-sample settings (provenance)
├── mk/                           # stage flag files (*.done) that drive make
├── tsv/
│   ├── in_list.tsv               # sample → transcript TSV (FICTURE input list)
│   ├── in_{boundaries,xy,clust}.tsv   # cell-decode input lists (when applicable)
│   └── <id>/transcripts.unsorted.tsv.gz
├── fic/                          # FICTURE results; per-sample under fic/samples/<id>/
└── cartl/samples/<id>/           # packaged PMTiles + catalog.yaml (per sample)
```

For a `--out-root` batch, this layout is created under `<out_root>/<id>/`, and the master Makefile and resolved config are written at `<out_root>/`.

!!! tip "Inspect before running"
    Run with `--dry-run` first and read `run_together.mk` and `run_together.resolved.json`. The Makefile shows the exact sub-commands and their dependencies; the resolved JSON shows the settings each sample was expanded to after all merges.

See the individual module references for the format of each output ([`sge_convert`](./sge_convert.md#output), [`run_ficture2_multi`](./run_ficture2_multi.md#output), [`run_cartload2`](./run_cartload2.md#output)).
