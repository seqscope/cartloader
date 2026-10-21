# Supported Platforms & Inputs

`run_together`'s one-command convenience comes from **platform profiles** that read each platform's standard output layout automatically. This page is the **status hub**; each built-in profile has its own page with the exact expected inputs and example commands.

For the profile mechanics (input modes, the canonical config, layer merging), see [Specifying Inputs](./run_together_inputs.md).

---
## Profile status

| `--platform` | Profile | Ingest | Details |
|--------------|:-------:|:------:|---------|
| `10x_xenium` | ✅ built-in | `sge_convert` (`10x_xenium`) | [Xenium page](./platforms/xenium.md) |
| `10x_visium_hd` | ✅ built-in | `sge_convert` (`10x_visium_hd`) | [Visium HD page](./platforms/visium_hd.md) |
| `cosmx_smi` | ✅ built-in | `reformat_cosmx` | [CosMx SMI page](./platforms/cosmx_smi.md) |
| `merfish` (Vizgen MERSCOPE) | ✅ built-in | `sge_convert` (`vizgen_merscope`) | [MERSCOPE page](./platforms/merscope.md) |
| `illumina` (Illumina StrataMap) | ✅ built-in | `sge_convert` (`illumina`) | [Illumina StrataMap page](./platforms/illumina.md) |
| `stereoseq` (BGI Stereo-seq) | ✅ built-in | SAW `gef2gem` + `sge_convert` (`bgi_stereoseq`) | [Stereo-seq page](./platforms/stereoseq.md) |
| `seqscope` (Seq-Scope) | ✅ built-in | `sge_convert` (`seqscope`) | [Seq-Scope page](./platforms/seqscope.md) |
| `generic` (any source) | ✅ built-in | an already-ingested transcript TSV, or `sge_convert` (`generic`) | [Generic page](./platforms/generic.md) |
| Pixel-seq, Nova-ST | via a custom profile | `sge_convert` (`pixel_seq`, `nova_st`) | [below](#custom-profiles) |

- **built-in** — `run_together --platform <name>` works with no extra configuration. The profiles are `assets/run_together_profiles/<name>.json`.
- **via a custom profile** — `sge_convert` already has an ingest preset; there is no dedicated `run_together` profile yet.

---
## Built-in platforms

Each page documents the expected input files/columns, the FICTURE/packaging defaults, and example commands for the input modes (single-sample CLI, sample sheet, config JSON):

- [**10x Xenium**](./platforms/xenium.md) — Xenium Ranger directory; morphology images auto-detected; `xeniumranger` clusters.
- [**10x Visium HD**](./platforms/visium_hd.md) — Space Ranger `outs/`; square-bin and segmented-cell layers; a `cartloader` cell analysis on the segmented cells; per-sample H&E.
- [**CosMx SMI**](./platforms/cosmx_smi.md) — AtoMx flat-file export via `reformat_cosmx`; glob-based file matching with overrides.
- [**MERSCOPE / MERFISH**](./platforms/merscope.md) — individual files or a standard export dir; cell analysis from boundaries, a cell×gene matrix, or an existing `cell_id` (incl. mixed joint runs).
- [**Illumina StrataMap**](./platforms/illumina.md) — MEX input with spatial coordinates embedded in the barcode; optional cell boundaries enable cell-level analysis.
- [**BGI Stereo-seq**](./platforms/stereoseq.md) — prefix-addressed binary GEFs expanded by SAW (`--saw`); cell-bin segmentation clustered from its own pixel TSV; registered histology at 0.5 µm/pixel.
- [**Seq-Scope**](./platforms/seqscope.md) — `--in-dir` *is* the MEX directory (or a raw per-molecule TSV via `--raw-transcript`); nanometer coordinates; pixel-level only (no segmentation); pre-georeferenced H&E.
- [**Generic**](./platforms/generic.md) — a transcript TSV from any source (already in CartLoader's layout, or a CSV/TSV converted by `sge_convert --platform generic`); optional cell analysis from a cell-id column or a MEX.

---
## Custom profiles

`--platform-json <file>` **deep-merges** a JSON file over the built-in profile named by `--platform`: nested blocks merge key by key, while lists (e.g. `ingest.extra_flags`, `ingest.autodetect`, `images`) are replaced whole. `--platform` must still name a built-in profile, so a custom profile is always a modification of one — start from the closest match, or from `generic` for a plain transcript table. The file is not checked for unknown keys.

For example, a Pixel-seq transcript table through `sge_convert`'s `pixel_seq` preset:

```json
{ "ingest": { "sge_platform": "pixel_seq" } }
```
```bash
cartloader run_together --platform generic --platform-json pixelseq.json \
    --raw-transcript /data/pixelseq/transcripts.tsv.gz --out-dir OUT
```

`ingest.sge_platform` selects the `sge_convert` preset (column names, delimiter); add `--units-per-um` if the coordinates are not in microns. Any profile key can be set this way; a `--config` can carry the same blocks for a single run (see [Specifying Inputs → Mode 3](./run_together_inputs.md#mode-3-full-config-json)), with strict key validation.

Without `run_together`, the modules can also be run directly: `sge_convert` → [`run_ficture2`](./run_ficture2.md) → [`run_cartload2`](./run_cartload2.md) (see the platform tutorials and [`sge_convert`](./sge_convert.md) for the per-platform ingest presets).

---
## See also

- [`run_together` Overview](./run_together.md) · [Specifying Inputs](./run_together_inputs.md) · [Image Modalities](./run_together_images.md)
- Guides: [single sample](../vignettes/pipelines/run_together.md) · [multiple samples](../vignettes/pipelines/run_together_multi.md)
- [`sge_convert` reference](./sge_convert.md) — per-platform ingest details.
