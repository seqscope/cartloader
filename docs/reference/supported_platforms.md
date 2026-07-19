# Supported Platforms & Inputs

`run_together`'s one-command convenience comes from **platform profiles** that read each platform's standard output layout automatically. This page is the **status hub**; each built-in platform has its own page with the exact expected inputs and example commands for all three input modes.

For the profile mechanics (input modes, the canonical config, layer merging), see [Specifying Inputs](./run_together_inputs.md).

---
## Profile status

| `--platform` | Profile | Ingest | Details |
|--------------|:-------:|:------:|---------|
| `10x_xenium` | ✅ built-in | `sge_convert` | [Xenium page](./platforms/xenium.md) |
| `10x_visium_hd` | ✅ built-in | `sge_convert` | [Visium HD page](./platforms/visium_hd.md) |
| `cosmx_smi` | ✅ built-in | `reformat_cosmx` | [CosMx SMI page](./platforms/cosmx_smi.md) |
| `merfish` (Vizgen MERSCOPE) | ✅ built-in | `sge_convert` (`vizgen_merscope`) | [MERSCOPE page](./platforms/merscope.md) |
| `illumina` (Illumina StrataMap) | ✅ built-in | `sge_convert` (`illumina`) | [Illumina StrataMap page](./platforms/illumina.md) |
| `bgi_stereoseq` (BGI Stereo-seq) | ✅ built-in | SAW `gef2gem` + `sge_convert` (`bgi_stereoseq`) | [Stereo-seq page](./platforms/stereoseq.md) |
| `seqscope` | ⏳ planned | `sge_convert` | [below](#platforms-without-a-built-in-profile-yet) |
| `pixel_seq`, `nova_st` | ⏳ planned | `sge_convert` | [below](#platforms-without-a-built-in-profile-yet) |
| `generic` (custom CSV/TSV) | via custom profile | `sge_convert` | [below](#platforms-without-a-built-in-profile-yet) |

- **built-in** — `run_together --platform <name>` works with no extra configuration.
- **planned** — `sge_convert` already supports ingest; a `run_together` profile is not yet shipped.

---
## Built-in platforms

Each page documents the expected input files/columns, the FICTURE/packaging defaults, and example commands for the three input modes (single-sample CLI, sample sheet, config JSON):

- [**10x Xenium**](./platforms/xenium.md) — Xenium Ranger directory; morphology images auto-detected; `xeniumranger` clusters.
- [**10x Visium HD**](./platforms/visium_hd.md) — Space Ranger `outs/`; automatic ×2 scaling; square + segmented-cell layers; per-sample H&E.
- [**CosMx SMI**](./platforms/cosmx_smi.md) — AtoMx flat-file export via `reformat_cosmx`; glob-based file matching with overrides.
- [**MERSCOPE / MERFISH**](./platforms/merscope.md) — individual files or a standard export dir; cell analysis from boundaries, a cell×gene matrix, or an existing `cell_id` (incl. mixed joint runs).
- [**Illumina StrataMap**](./platforms/illumina.md) — MEX input with spatial coordinates embedded in the barcode; optional cell boundaries enable cell-level analysis.
- [**BGI Stereo-seq**](./platforms/stereoseq.md) — prefix-addressed binary GEFs expanded by SAW (`--saw`); cell-bin segmentation clustered from its own pixel TSV; registered histology at 0.5 µm/pixel.

---
## Platforms without a built-in profile yet

`sge_convert` converts every platform in the status table to the unified transcript TSV that FICTURE and packaging consume — the per-platform column names, delimiters, and scaling are already encoded there (see [`sge_convert`](./sge_convert.md)). Two ways to run them end-to-end today:

1. **Modules directly** — `sge_convert` → [`run_ficture2`](./run_ficture2.md) → [`run_cartload2`](./run_cartload2.md). See the platform starter tutorials for worked examples.
2. **A custom `run_together` profile** — copy a built-in profile (`assets/run_together_profiles/10x_xenium.json`), set `ingest.sge_platform` to the target platform and adjust the input paths / exclusion regex, then pass it with `--platform-json`. See [Specifying Inputs → Mode 3](./run_together_inputs.md#mode-3-full-config-json).

---
## See also

- [`run_together` Overview](./run_together.md) · [Specifying Inputs](./run_together_inputs.md) · [Image Modalities](./run_together_images.md)
- [Single-sample tutorial](../vignettes/pipelines/run_together.md) · [multi-sample tutorial](../vignettes/pipelines/run_together_multi.md)
- [`sge_convert` reference](./sge_convert.md) — per-platform ingest details.
