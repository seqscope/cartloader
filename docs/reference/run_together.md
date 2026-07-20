# Unified Pipeline Orchestrator (`run_together`)

## Overview

`run_together` runs a complete CartoScope pipeline — **ingest → FICTURE → cell decode → asset packaging → image import → (optional) publish** — across many spatial platforms, for one or many samples, from a single command. It assembles the whole run as **one Makefile** (`run_together.mk`) and executes it, so the pipeline is **resumable**, **parallel**, and trains **one joint FICTURE model** across a multi-sample run.

The design has one guiding principle:

> **Convenience for the default case; full control in JSON.**

A standard run needs only a few flags. Everything a run *can* express lives in a single canonical, list-based configuration that three layers assemble:

1. **Profile** — per-platform defaults (auto-detected inputs, default analyses, image conventions).
2. **Tier-1 CLI** — the common knobs (`--in-dir`/`--samples`, `--width`/`--n-factor` or `--project-models`, `--out-dir`).
3. **Tier-2 JSON** (`--config`) — augments or fully specifies anything.

The layers **compose**: set the base run on the CLI and add only the extra analyses/images in JSON.

### Read next

| To… | See |
|-----|-----|
| Specify inputs (single CLI, sample sheet, or config JSON) | [**Specifying Inputs**](./run_together_inputs.md) |
| Attach DAPI / H&E / protein images | [**Image Modalities**](./run_together_images.md) |
| See exact inputs + examples for your platform | [**Platforms**](./supported_platforms.md): [Xenium](./platforms/xenium.md) · [Visium HD](./platforms/visium_hd.md) · [CosMx SMI](./platforms/cosmx_smi.md) · [MERSCOPE](./platforms/merscope.md) |

A minimal run:

```bash
cartloader run_together --platform 10x_xenium --in-dir IN --out-dir OUT \
    --width 12 --n-factor 24
```

---
## Requirements

- `make` on the `PATH`.
- Tools used by the stages: `spatula`, `punkst` (FICTURE2), `pigz`, `sort`, `python3`, `go-pmtiles`, `gdal`, `tippecanoe`, `parquet-tools`, `jq` (Visium HD H&E), and `aws` (publish). `run_together` delegates to the CartLoader modules, which find `spatula`/`punkst` as built repo submodules or on `PATH`.

---
## Stages, selection & resume

The pipeline stages are `ingest, ficture, cells, cartload, images, anno, upload`.

- **Resume:** re-run the same command (or `make -f OUT/run_together.mk -j N`); completed stages are skipped via flag files.
- `--only ingest,ficture` / `--skip images` — run a subset; excluded upstream stages are assumed done (prereqs are pruned so `make` won't error).
- `--restart` — rebuild everything (`make -B`).
- `--dry-run` — write the Makefile and print commands (`make -n`) without executing.

**Packaging bundles everything produced** — every FICTURE pixel decode plus every cell analysis that ran.

---
## Publishing

Publishing is **opt-in** and **entirely CLI-driven** (no config block). Two independent actions — enable either or both:

- **`--anno`** — AI-annotate each packaged sample directory (via [`anno_cartload_folder`](./anno_cartload_folder.md)). **Requires `--tissue` and `--organism`** (no defaults); `--anno-api-type` (default `umgpt`), `--anno-model` (default `claude-opus-4-7`), and `--anno-threads` (default `10`) are overridable. For a joint run the shared factors are annotated **once** at the `cartl/` root and reused into every sample.
- **`--s3-upload`** — upload each self-contained sample directory to S3.

When both run, `upload` waits on `anno` (which edits `catalog.yaml`).

**S3 destination:** `<s3-prefix>/batch=<batch>/<collection>/<dir-id>/`, where `<dir-id>` is the sample's output directory name (`<sample_id>` single, `<multi_id>-<sample_id>` joint). For a **joint run**, `multi-catalog.yaml` and the shared factor files are additionally uploaded to the **parent** `<s3-prefix>/batch=<batch>/<collection>/` so relative pointers resolve.

- `--s3-prefix` — default `s3://cartostore/data` · `--batch` — default current `YYYY_MM` · `--collection` — default out-dir basename
- `--aws-profile` (default `cartostore`) · `--aws` (binary path) · `--s3-jobs` (parallel copies, default `4`)

```bash
cartloader run_together --platform 10x_xenium --samples samples.tsv --out-dir OUT \
    --width 12 --n-factor 24 \
    --anno --tissue "Kidney" --organism human \
    --s3-upload --collection my-collection
```

---
## How samples are packaged

- **Single sample** → one `run_cartload2` call, output at `cartl/<id>/`.
- **Joint multi-sample run** (samples sharing one `--out-dir`) → a single [`run_cartload2_multi`](./run_cartload2_multi.md) call that packages every sample in parallel. Each sample is a **self-contained** `cartl/<multi_id>-<sample_id>/` directory, and a **`cartl/multi-catalog.yaml`** links every per-sample `catalog.yaml` and copies the shared factor files (`post`/`rgb`/`de`/`info`/`umap`) into the `cartl/` root — so `cartl/` uploads to S3 as one deployable unit.

This mirrors the FICTURE manifests: `run_ficture2_multi` writes a shared [`ficture.multi.params.json`](./run_ficture2_multi.md#multi-sample-manifest) that `run_cartload2_multi` reads to discover samples and shared assets.

---
## Command-line parameters

**Run:** `--dry-run`, `--restart`, `-j/--n-jobs`, `--threads`, `--makefn`, `--only`, `--skip`.

**Input/output:** `--platform`, `--in-dir`, `--samples`, `--out-dir`, `--out-root`, `--id`, `--config`, `--platform-json` (external profile override). Single-sample file inputs and column-name overrides are documented in [Specifying Inputs](./run_together_inputs.md).

**FICTURE mode:** `--width`, `--n-factor` (de-novo); `--project-models` (projection-only); `--no-ficture` (tiling only — package with no factor layers). See [Specifying Inputs → FICTURE mode](./run_together_inputs.md#ficture-mode-de-novo-vs-projection).

**Common decode overrides:** `--exclude-feature-regex`, `--min-ct-per-unit-hexagon` (default `50`), `--always-single-molecule` / `--never-single-molecule` (default: single-molecule ON for pixel FICTURE, OFF for cell decode). An explicit CLI flag wins over a `--config`/profile value, which wins over the built-in default.

**Publish:** `--anno`, `--s3-upload`, `--tissue`, `--organism`, `--anno-api-type`, `--anno-model`, `--anno-threads`, `--collection`, `--batch`, `--s3-prefix`, `--aws-profile`, `--aws`, `--s3-jobs` (see [Publishing](#publishing)).

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

For a joint run, upload the whole `cartl/` directory. For `--out-root` batches (independent per-sample models), this layout is created under `<out_root>/<id>/`, and the master Makefile is written at `<out_root>/`.

---
## See also

- [Specifying Inputs](./run_together_inputs.md) · [Image Modalities](./run_together_images.md) · [Supported Platforms](./supported_platforms.md)
- Tutorials: [single-sample](../vignettes/pipelines/run_together.md) · [multi-sample](../vignettes/pipelines/run_together_multi.md)
