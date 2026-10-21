# Annotate Factors with an LLM (`annotate_factors_with_llm`)

!!! tip "Usually run for you"
    [`run_together --anno-deep`](./run_together.md#publishing) runs this command through [`anno_cartload_folder --deep`](./anno_cartload_folder.md), once per factor model, and records the result in the catalog. This page is for running the annotator directly, step by step.

## Overview

`annotate_factors_with_llm` labels the latent factors of a cartloader run (cell types or transcriptional programs) with a generative LLM that sees **all factors at once**, so it can contrast factors and answer `Unresolved` when the evidence is weak. For each factor it returns a best-guess label, kind (cell type or program), compartment, confidence, key genes, rationale, and alternative interpretations when the markers support more than one reading. It writes an alias TSV in the same format as `annotate_bulk_de_with_ai` (a drop-in for the catalog's `alias_ai` key), a per-factor annotation table, and a self-contained interactive HTML report.

Evidence per factor:

- marker genes: the top `--top-n` genes by chi-square and by fold enrichment, interleaved;
- abundance (share of counts), DE support (number of enriched genes), and per-gene specificity (needs the model matrix);
- multi-sample: which samples and tissues hold the factor, and its top genes within each sample.

All factors are sent in one call by default; a call whose answer exceeds the output budget is split in halves automatically.

---
## Requirements

- Packaged factor files, addressed by a prefix (see [Inputs](#inputs)): a model matrix (`-model.tsv`) and/or a DE table (`-bulk-de.tsv` or `-de.tsv`; computed from the model when absent).
- An API key in the environment for `--api-type`: `ANTHROPIC_API_KEY` (or `ANTHROPIC_AUTH_TOKEN`) for `claude`, `OPENAI_API_KEY` for `openai`, `GEMINI_API_KEY` for `google`, `UMGPT_API_KEY` for `umgpt`. Not needed with `--dry-run` or `--report-only`.
- For `--api-type claude`: the `anthropic` Python package, version 0.77.0 or later.

---
## Inputs

Point `--prefix` at a factor prefix of a packaged directory:

| Run | Example prefix | Files read |
|---|---|---|
| Single sample | `cartl/<sample>/t18-f192` | `<prefix>-bulk-de.tsv`, `<prefix>-model.tsv` |
| Multi-sample root | `cartl/t12-f192` | `<prefix>-de.tsv`, `<prefix>-model.tsv`, and each sample's `<multi_id>-<sample>/t12-f192-pixel-pseudobulk.tsv.gz` |

`<prefix>-rgb.tsv`, when present, gives each factor its CartoScope color in the report. For the factor model, `-model.tsv`, `-model.tsv.gz`, or `-pseudobulk.tsv.gz` is used, whichever exists first.

Samples of a multi-sample run are found through `multi-catalog.yaml` next to the prefix or, without it, from sibling sub-directories holding a pseudobulk for the prefix. A run with a single sample is annotated as single-sample. The tissue of every sample is `--tissue`, or the `tissue` column of `--sample-sheet` (e.g., the `run_together` sample sheet, matched by `id`) for pan-tissue runs. A single-sample dataset requires `--tissue`.

---
## Example Usage

```bash
# Single sample
cartloader annotate_factors_with_llm \
  --prefix cartl/sample-rep1/t18-f192 \
  --organism mouse --tissue "whole mouse pup" \
  --api-type claude

# Multi-sample root with per-sample tissues
cartloader annotate_factors_with_llm \
  --prefix cartl/t12-f192 \
  --organism human --sample-sheet samples.tsv \
  --api-type claude

# Preview the prompts and an evidence-only report without calling the LLM
cartloader annotate_factors_with_llm \
  --prefix cartl/sample-rep1/t18-f192 \
  --organism mouse --tissue "whole mouse pup" \
  --api-type claude --dry-run
```

---
## Candidate labels

Candidate labels are optional suggestions sent with the prompt (or a closed list with `--closed`).

- `--candidates auto` generates one list per tissue with the same provider (reasoning effort `--candidate-effort`, default `low`), caches it under `--candidate-cache-dir`, and saves the merged list as `<stem>.candidates.json`. Uncached lists are generated with an LLM call even under `--dry-run`; `--report-only` reuses the saved `<stem>.candidates.json`.
- `--candidates FILE [FILE ...]` uses JSON lists made by `generate_candidate_labels_with_llm`; multiple files are merged by label.

`generate_candidate_labels_with_llm` makes one LLM call per tissue context (no marker genes are sent), so a list can be reused for every dataset from that tissue:

```bash
cartloader generate_candidate_labels_with_llm \
  --organism human --tissue "lung cancer" "lymph node" \
  --api-type claude --out candidates.json
```

Its options: `--organism`, `--tissue` (one or more), `--out`, `--api-type` (all required); `--mode` (default `either`), `--cache-dir` (default `./jevanno_cache`), `--no-generate` (use cached lists only), `--min-candidates` / `--max-candidates` (default `80` / `200`), `--model-name`, `--effort` (default `low`), `--max-output-tokens` (default `64000`), `--threads` (default `8`), plus `--no-fallbacks`, `--api-base-url`, `--request-timeout`, `--verbose`. Lists are cached per model and effort, and a cached list made at the same or higher effort is reused.

---
## Parameters

### Input/Output Parameters

- `--prefix` (str): Factor prefix (see [Inputs](#inputs)).
- `--model` (str): Factor model matrix TSV, gene × factor (default: from `--prefix`).
- `--de` (str): Bulk DE TSV (default: from `--prefix`; computed from the model when absent).
- `--rgb` (str): Factor color table (default: `<prefix>-rgb.tsv` if present).
- `--out` (str): Output alias TSV (default: `<prefix>-alias-ai.tsv`); companion files share its stem. Required without `--prefix`.
- `--organism` (str, required): e.g., `human`, `mouse`.
- `--tissue` (str): Tissue context for all samples (unless the sample sheet gives one per sample).
- `--sample-sheet` (alias `--samples`) (str): Multi-sample TSV with `id` and optional `tissue` (and `pseudobulk`) columns.
- `--candidates` (list of str): Candidate label JSON file(s), or `auto`.
- `--api-type` (str, required): `claude`, `openai`, `umgpt`, or `google`.

??? note "Auxiliary Parameters"

    **Annotation**:

    - `--mode` (str): `cell_type`, `program`, or `either` (default: `either`).
    - `--closed` (flag): Restrict labels to the candidates (requires `--candidates`).
    - `--single-sample` (flag): Treat the prefix as single-sample even if sample directories are found (cannot be combined with `--sample-sheet`).
    - `--multi-catalog` (str): Multi-sample catalog next to the prefix (default: `multi-catalog.yaml`).
    - `--candidate-cache-dir` (str): Cache for `--candidates auto` (default: `~/.cache/cartloader/llm_candidates`).
    - `--candidate-effort` (str): Effort for `--candidates auto` (default: `low`).

    **Evidence**:

    - `--primary-rank` / `--secondary-rank` (str): DE ranking columns (default: `Chi2` / `FoldChange`).
    - `--top-n` (int): Top genes per ranking column (default: `10`).
    - `--min-marker-gene-count` (int): Factors with fewer marker genes are left `Unresolved` without an API call (default: `5`).
    - `--no-specificity` (flag): Omit per-gene specificity.
    - `--no-sample-markers` (flag): Multi-sample: omit per-sample marker genes.

    **LLM**:

    - `--model-name` (str): LLM model (defaults: `claude-opus-5-5` for `claude`, `gpt-5.4-mini` for `openai`, `gpt-5-mini` for `umgpt`, `gemini-3.5-flash` for `google`).
    - `--effort` (str): Reasoning effort: `low`, `medium`, `high`, `xhigh`, `max` (default: `high`).
    - `--max-factors-per-call` (int): Factors per call; `0` sends all factors in one call (default: `0`).
    - `--max-output-tokens` (int): Output budget per call (defaults: `128000` for `claude`/`openai`, `64000` for `umgpt`/`google`).
    - `--threads` (int): Parallel LLM calls when there are several (default: `4`).
    - `--request-timeout` (int): Seconds (default: `1800`); `--max-retries` (int): Retries (default: `3`).
    - `--api-base-url` (str): Base URL for `openai`/`umgpt` (e.g., a proxy).
    - `--no-fallbacks` (flag): `claude`: disable server-side refusal fallbacks.

    **Run and Report**:

    - `--dry-run` (flag): Write the prompts to `<stem>.llm/` and an evidence-only report, without calling the LLM (no TSVs are written).
    - `--report-only` (flag): Rebuild the report and TSVs from responses saved in `<stem>.llm/` without calling the LLM.
    - `--no-reuse` (flag): Do not reuse responses saved in `<stem>.llm/` (they are reused by default).
    - `--html` (str): Report path (default: `<stem>.html`); `--no-html` (flag) skips the report.
    - `--dotplot-genes` (str): Genes embedded in the report's dotplots: `markers` (default) or `all` (larger file).
    - `--title` (str): Report title (default: dataset and prefix name).
    - `--verbose` (flag): Verbose logging.

---
## Output

With the default stem `<prefix>-alias-ai`:

- `<stem>.tsv`: Alias TSV (`index`, `alias`); duplicate aliases get `_1`, `_2`, … suffixes and unresolved factors are `Unresolved`.
- `<stem>.annotations.tsv`: One row per factor with `status` (`assigned`, `ambiguous`, `unfitted`), alias, label, kind, compartment, confidence, key genes, alternatives, rationale, and evidence columns.
- `<stem>.html`: Self-contained report with an overview, the factor table, interactive dotplots, and run details.
- `<stem>.llm/`: Requests and responses of every LLM call (reused on re-runs unless `--no-reuse`).
- `<stem>.candidates.json`: The merged candidate list (with `--candidates auto`).
