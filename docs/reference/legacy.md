# Legacy Orchestrators

!!! warning "Superseded by `run_together`"
    `run_xenium` and `run_visiumhd` are the original single-sample, single-platform orchestrators. They are superseded by [`run_together`](./run_together.md), kept for reference only, and not part of the recommended workflow. New runs should use `run_together`, or the [individual commands](./intro.md#individual-commands) for step-by-step control.

## Overview

- [`run_xenium`](./run_xenium.md) runs one 10x Xenium sample through asset detection ([`load_xenium_ranger`](./xenium_addon.md)), SGE conversion, FICTURE, cell and image import, packaging, and upload.
- [`run_visiumhd`](./run_visiumhd.md) does the same for one 10x Visium HD sample ([`load_space_ranger`](./visiumhd_addon.md)), plus square-bin import.

Each action is a flag, and each flag calls one of the individual commands.

---
## What `run_together` covers

[`run_together`](./run_together.md) runs the same steps as one resumable Makefile, for one sample or many (with a joint FICTURE model), on every [supported platform](./supported_platforms.md), including [Xenium](./platforms/xenium.md) and [Visium HD](./platforms/visium_hd.md).

| Legacy action | `run_together` equivalent |
|---|---|
| `--load-xenium-ranger` / `--load-space-ranger` | The platform profile detects inputs under `--in-dir` |
| `--sge-convert` | `ingest` stage |
| `--run-ficture2` | `ficture` stage |
| `--import-cells` | `cells` stage (and per-platform cell import) |
| `--import-squares` (Visium HD) | Square-bin import in the Visium HD profile |
| `--import-images` | `images` stage ([Image Modalities](./run_together_images.md)) |
| `--run-cartload2` | `cartload` stage |
| `--upload-aws` | `--s3-upload` ([Publishing](./run_together.md#publishing)) |
| `--upload-zenodo` | None; run [`upload_zenodo`](./upload_zenodo.md) separately |

---
## Remaining gaps

Functions the legacy orchestrators have that `run_together` does not:

- **Zenodo upload** (`--upload-zenodo`). Run [`upload_zenodo`](./upload_zenodo.md) on the packaged directory instead.
- **Asset manifests and archive extraction.** `load_xenium_ranger` and `load_space_ranger` write a JSON manifest of detected files and extract `.tar.gz` archives (Xenium `analysis.tar.gz` and `cell_feature_matrix.tar.gz`; Space Ranger `*_binned_outputs.tar.gz` and `*_segmented_outputs.tar.gz`). `run_together` expects an unpacked `--in-dir`; the loaders can still be run on their own to unpack archives (see [Xenium Add-ons](./xenium_addon.md) and [Visium HD Add-ons](./visiumhd_addon.md)).
- **Catalog title and description** (`--title`, `--desc`, forwarded to `run_cartload2`). `run_together` sets only the catalog `id`, so the title defaults to the id.
- **Arbitrary S3 destination** (`--s3-dir`, via [`upload_aws`](./upload_aws.md)). `run_together` uploads to `<s3-prefix>/batch=<YYYY_MM>/<collection>/<dir>`.

---
## Known issues

The following problems were confirmed in the current code. In each case, running the individual command directly works around it.

- **`--run-ficture2` with `--run-cartload2` fails.** The packaging step waits for decode outputs named `t<w>_f<k>_p<w>_a6.done`, `..._summary.done`, and `.png` (`define_ficture_flags` in `cartloader/utils/pipeline_helper.py`), but `run_ficture2` now names them `t<w>_f<k>_pixel.done`, `t<w>_f<k>_pixel_summary.done`, and `t<w>_f<k>_pixel.png`, so the run stops with `Found missing prerequisite`. (`--dry-run` does not check prerequisites, so it does not show the problem.) Workaround: after FICTURE finishes, run `cartloader run_cartload2 --fic-dir <out-dir>/ficture2 --out-dir <out-dir>/cartload2 --id <id>` (add `--cell-assets` / `--background-assets` as needed).
- **`--import-ext-ficture2` fails.** It passes `--ext-fic-dir` to `run_cartload2`, which does not accept that option. Workaround: run `run_cartload2 --fic-dir <external FICTURE dir>`.
- **`--run-cartload2` without FICTURE fails.** It packages through `run_cartload2 --sge-dir`, which needs a tiled transcript that `sge_convert` does not produce (see [`run_cartload2`](./run_cartload2.md#requirements)). Workaround: run FICTURE first and package with `--fic-dir`, or use `run_together`.
- **`run_xenium --import-cells` in manual mode (`--csv-*`) fails.** `import_xenium_cell` regenerates pseudobulk and DE by default and needs a cell-feature MEX or a transcript TSV, which `run_xenium` does not pass in manual mode. Auto-detect mode works when the manifest includes `CELL_FEATURE_MEX`. Workaround: run [`import_xenium_cell`](./import_cell.md) directly with `--mex-dir cell_feature_matrix`.
- **`run_xenium --import-images` with the default image IDs fails** unless all five default images (`DAPI_OME`, `BOUNDARY_OME`, `INTERIOR_RNA_OME`, `INTERIOR_PROTEIN_OME`, `DAPI_MIP_OME`) are in the manifest. Pass `--image-ids` for the images you have, or `--all-images`.
- **`run_visiumhd --import-squares` with `--square-input`** does not pass the bin-level MEX directory, so [`import_visiumhd_square`](./import_square.md) looks for `<space-ranger-dir>/filtered_feature_bc_matrix`; it also requires the UMAP CSV. Use auto-detect mode, or run `import_visiumhd_square` directly with `--in-dir <space-ranger-dir>/binned_outputs/square_008um`.
- **`run_visiumhd --use-parquet-tools --parquet-tools <path>`** passes `--parquet-tools` without its value to `import_visiumhd_square`, which then fails to parse its arguments.
