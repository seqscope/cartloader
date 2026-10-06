# AI Annotation (`anno_cartload_folder`)

## Overview

`anno_cartload_folder` AI-annotates every factor of a packaged CartoScope dataset. For each factor listed in the catalog it (a) calls [`annotate_bulk_de_with_ai`](#) on the factor's DE table to produce a per-factor **alias TSV** (`<factor_id>-alias-ai.tsv`), and (b) records that alias file in the catalog under the `alias_ai` key. It operates either on a **local `--cartl-dir`** or an already-uploaded **`--s3-dir`** (downloading, annotating, and re-uploading).

The [`run_together`](./run_together.md) `--anno` stage invokes this script; you can also run it standalone on a packaged directory.

---
## Modes

=== "Single sample (default)"

    Annotates one packaged directory. The catalog is either a per-sample `catalog.yaml` (`assets.factors: [...]` list) or any top-level `factors: {…}` map.

    ```bash
    cartloader anno_cartload_folder \
      --cartl-dir cartl/mydataset \
      --tissue "Kidney" --organism human
    ```

=== "Multi-sample (`--multi-sample`)"

    Annotates the **shared** factors of a joint run **once** at the `cartl/` root (using the `multi-catalog.yaml` written by [`run_cartload2_multi`](./run_cartload2_multi.md)), then **propagates** each shared alias TSV into every sample sub-folder listed in `samples:`. Sub-folders reuse the shared results (via `--reuse-results-from`), so a factor is annotated at most once for the whole cohort and every sample's `catalog.yaml` ends up pointing at a matching alias file.

    ```bash
    cartloader anno_cartload_folder \
      --cartl-dir cartl \
      --multi-sample \
      --tissue "Kidney" --organism human
    ```

    Factor IDs are normalized (`t12_f48` → `t12-f48`) when naming alias files, so a shared alias matches the per-sample factor id in each `<multi_id>-<sample_id>/catalog.yaml`.

=== "S3 in-place (`--s3-dir`)"

    Downloads the catalog (and the per-factor DE files that need annotation), annotates locally, and uploads the alias TSVs and the updated catalog back to S3. Uses `boto3`; upload uses `--profile`.

    ```bash
    cartloader anno_cartload_folder \
      --s3-dir s3://cartostore/data/batch=2026_07/mycollection/mydataset \
      --tissue "Kidney" --organism human --profile cartostore
    ```

---
## Reusing prior annotations

Pass `--reuse-results-from <dir>` to skip the AI call and copy alias TSVs from `<dir>` instead. Useful for re-runs or for hand-annotating one sample and propagating its results to the rest of a cohort. In `--multi-sample` mode this is applied only to the shared (root) pass; the sub-folder propagation always reuses the shared results just written.

---
## Parameters

**Input/output (choose one):**

- `--cartl-dir` — local packaged directory
- `--s3-dir` — `s3://…` prefix of an uploaded dataset

**Required:**

- `--tissue` — tissue name (passed to the annotation prompt)
- `--organism` — organism/species

**Common:**

- `--multi-sample` — treat `--cartl-dir` as a joint-run root; use `multi-catalog.yaml` and propagate
- `--multi-catalog` — multi-catalog filename (default `multi-catalog.yaml`)
- `--reuse-results-from` — copy alias TSVs from this dir instead of annotating
- `--api-type` (default `umgpt`), `--model` (default `claude-opus-4-7`), `--threads` (default `1`)
- `--catalog` — per-sample catalog filename (default `catalog.yaml`)

**Alias / catalog keys:**

- `--alias-suffix` (default `-alias-ai.tsv`) — suffix for alias file names
- `--backup-suffix` (default `.bak`) — appended to the pre-update catalog copy
- `--yaml-key-store` (default `alias_ai`) — factor key that records the alias file
- `--yaml-key-skip` (default `alias alias_ai`) — factor keys that mark an already-annotated factor (skip if present)

**S3 (with `--s3-dir`):**

- `--profile` (default `cartostore`), `--aws` (default `aws`)
- `--skip-upload` — annotate locally in the temp dir, do not upload back
- `--tmp-dir` — where to download / stage files

---
## Output

Under the target directory:

- `<factor_id>-alias-ai.tsv` for every factor that had a `de` entry and no prior alias.
- The catalog rewritten in place, with `alias_ai: <factor_id>-alias-ai.tsv` added to each annotated factor.
- A backup of the pre-update catalog at `<catalog><backup-suffix>` (default `.bak`).

For `--multi-sample`: alias files are written both at the `cartl/` root and inside every sample sub-folder listed in `multi-catalog.yaml`'s `samples:`.
