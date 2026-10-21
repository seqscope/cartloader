# AI Annotation (`anno_cartload_folder`)

!!! tip "Usually run for you"
    [`run_together`](./run_together.md#publishing) calls `anno_cartload_folder` in its `anno` stage: `--anno` uses the per-factor annotator and `--anno-deep` adds `--deep`. This page is for annotating a packaged directory directly, step by step.

## Overview

`anno_cartload_folder` AI-annotates every factor of a packaged CartoScope dataset. For each factor listed in the catalog it writes an **alias TSV** (`<factor-id>-alias-ai.tsv`, factor index → label) and records that file in the catalog under the `alias_ai` key. Two annotators are available:

- **Default**: `annotate_bulk_de_with_ai` labels each factor from its DE table, one LLM call per factor.
- **`--deep`**: [`annotate_factors_with_llm`](./annotate_factors_with_llm.md) labels all factors of a model together in one deep-reasoning call, with alternative interpretations, and also writes an annotation table and an interactive HTML report.

It operates on a **local `--cartl-dir`** or an already-uploaded **`--s3-dir`** (download, annotate, upload back).

---
## Requirements

- A packaged directory with a `catalog.yaml` (from [`run_cartload2`](./run_cartload2.md)) or, for `--multi-sample`, a `run_cartload2_multi` output root with `multi-catalog.yaml`.
- An API key for the chosen `--api-type` in the environment: `ANTHROPIC_API_KEY` (`claude`), `OPENAI_API_KEY` (`openai`), `GEMINI_API_KEY` (`google`), or `UMGPT_API_KEY` (`umgpt`).
- For `--s3-dir`: the `boto3` package and AWS credentials.

---
## Modes

=== "Single sample (default)"

    Annotates one packaged directory, reading `assets.factors` from its `catalog.yaml`.

    ```bash
    cartloader anno_cartload_folder \
      --cartl-dir cartl/mydataset \
      --tissue "Kidney" --organism human
    ```

=== "Multi-sample (`--multi-sample`)"

    Annotates the **shared** factors of a joint run **once** at the `cartl/` root (from the `multi-catalog.yaml` written by [`run_cartload2_multi`](./run_cartload2_multi.md)), then copies each shared alias TSV into every sample sub-folder listed under `samples:` and records it in that sample's `catalog.yaml`. A sample-specific factor without a shared annotation is annotated in its own folder.

    ```bash
    cartloader anno_cartload_folder \
      --cartl-dir cartl \
      --multi-sample \
      --tissue "Kidney" --organism human
    ```

    Factor IDs are hyphenated (`t12_f48` → `t12-f48`) in alias file names, so a shared alias matches the per-sample factor id in each `<multi_id>-<sample_id>/catalog.yaml`.

=== "Deep annotation (`--deep`)"

    Runs `annotate_factors_with_llm` per factor model. In a multi-sample root, per-sample tissues can come from a sample sheet (`id` and `tissue` columns), which may replace `--tissue`.

    ```bash
    cartloader anno_cartload_folder \
      --cartl-dir cartl \
      --multi-sample --deep \
      --sample-sheet samples.tsv \
      --organism human --api-type claude --effort high
    ```

=== "S3 in-place (`--s3-dir`)"

    Downloads the catalog and the files needed for annotation, annotates them in a local temporary directory, and uploads the alias TSVs and the updated catalog back to S3. Downloads use the default `boto3` credentials; uploads use `--profile`. `--multi-sample` has no effect in this mode.

    ```bash
    cartloader anno_cartload_folder \
      --s3-dir s3://cartostore/data/batch=2026_07/mycollection/mydataset \
      --tissue "Kidney" --organism human --profile cartostore
    ```

---
## Skipping and reusing annotations

- Factors that already carry an `alias` or `alias_ai` key (`--yaml-key-skip`) are skipped, so reruns only annotate new factors. Factors without a `de` entry are skipped, except in `--deep` mode when they have a model (`model` or `post`) entry.
- Pass `--reuse-results-from <dir>` to copy alias TSVs from `<dir>` instead of calling the AI, e.g., to propagate one hand-checked sample to a cohort. In `--multi-sample` mode it applies only to the shared (root) pass; the sample propagation always reuses the shared results just written.

---
## Parameters

**Input (choose one):**

- `--cartl-dir` (str): Local packaged directory.
- `--s3-dir` (str): `s3://…` prefix of an uploaded dataset.

**Annotation context:**

- `--organism` (str, required): Organism/species.
- `--tissue` (str): Tissue name. Required, except with `--deep --multi-sample` when `--sample-sheet` gives samples a tissue (samples without one fall back to `--tissue`).

**Annotator:**

- `--api-type` (str): `claude`, `openai`, `google`, or `umgpt` (default: `claude`).
- `--model` (str): LLM model name (default: `claude-opus-4-8`; with `--deep`, the `annotate_factors_with_llm` default for the API type, e.g. `claude-opus-5-5` for `claude`).
- `--threads` (int): Annotation threads (default: `1`).
- `--deep` (flag): Use `annotate_factors_with_llm` instead of `annotate_bulk_de_with_ai`.
- `--effort` (str): `--deep` reasoning effort: `low`, `medium`, `high`, `xhigh`, or `max` (default: `high`).
- `--sample-sheet` (str): `--deep` with `--multi-sample`: TSV with `id` and `tissue` columns (blank, `-`, `.`, or `NA` mean unset).

**Multi-sample and reuse:**

- `--multi-sample` (flag): Treat `--cartl-dir` as a joint-run root; annotate `multi-catalog.yaml` and propagate.
- `--multi-catalog` (str): Multi-catalog file name (default: `multi-catalog.yaml`).
- `--catalog` (str): Per-sample catalog file name (default: `catalog.yaml`).
- `--reuse-results-from` (str): Copy alias TSVs from this directory instead of annotating.

**Alias and catalog keys:**

- `--alias-suffix` (str): Alias file suffix (default: `-alias-ai.tsv`).
- `--backup-suffix` (str): Suffix of the pre-update catalog copy (default: `.bak`).
- `--yaml-key-store` (str): Factor key that records the alias file (default: `alias_ai`).
- `--yaml-key-skip` (list of str): Factor keys that mark an annotated factor (default: `alias alias_ai`).

**S3 (with `--s3-dir`):**

- `--profile` (str): AWS profile for uploading results (default: `default`).
- `--skip-upload` (flag): Annotate locally only; do not upload back.
- `--tmp-dir` (str): Download and staging directory (default: a directory named after the last component of `--s3-dir`, under the current directory).

---
## Output

In each annotated directory:

- `<factor-id>-alias-ai.tsv`: Alias TSV (`index`, `alias`) for every annotated factor.
- With `--deep`, next to it: `<factor-id>-alias-ai.annotations.tsv` (per-factor label, confidence, alternatives, evidence), `<factor-id>-alias-ai.html` (report), and `<factor-id>-alias-ai.llm/` (saved LLM requests and responses).
- The catalog rewritten in place with `alias_ai: <factor-id>-alias-ai.tsv` on each annotated factor, plus a backup of the previous catalog at `<catalog><backup-suffix>` (written on every run).

For `--multi-sample`, alias files are written at the `cartl/` root and in every sample sub-folder listed under `samples:` in `multi-catalog.yaml`.
