# Feature Customization

## Overview

Features (typically genes) can be filtered at several points of the pipeline. Most commands take the same four flags — `--include-feature-list`, `--exclude-feature-list`, `--include-feature-regex`, `--exclude-feature-regex` — but what they affect depends on the command:

| Command | Filter flags | What is filtered |
|---------|--------------|------------------|
| [`sge_convert`](./sge_convert.md#feature-filtering-parameters) | lists and regexes | The transcript TSV at ingest, so everything downstream |
| [`reformat_cosmx`](./reformat_cosmx.md#feature-filtering) | lists and regexes | The CosMx transcript TSV at ingest |
| [`run_ficture2`](./run_ficture2.md#key-parameters) | regexes; `--in-feature-ficture` | The LDA model only (the tiled data keep every feature) |
| [`run_ficture2_multi`](./run_ficture2_multi.md#feature-filters) | lists and regexes | The LDA training/projection only, via [`feature_select`](#resolving-filters-into-one-list-feature_select) |
| [`run_ficture2_multi_cells`](./run_ficture2_multi_cells.md#feature-filters) | lists and regexes | The per-cell count matrices, so the whole cell analysis |
| [`feature_filtering`](#filtering-a-feature-table-feature_filtering) | lists, substrings, regexes, types | Rows of any CSV/TSV with a feature column (standalone) |
| [`feature_select`](#resolving-filters-into-one-list-feature_select) | lists and regexes | Writes one resolved feature list |

!!! tip "Usually run for you"
    [`run_together`](./run_together.md) passes its ingest filters to `sge_convert`/`reformat_cosmx` and its analysis filters to `run_ficture2_multi`, which calls `feature_select` (see [Feature filtering: two independent layers](./run_together_inputs.md#feature-filtering-two-independent-layers)). This page is for running the filters directly in step-by-step or custom workflows; `feature_filtering` is a standalone add-on that `run_together` does not call.

### Regex semantics differ between commands

- **Substring (search)**: `sge_convert`, `reformat_cosmx`, `run_ficture2_multi`, `run_ficture2_multi_cells`, `feature_filtering`, and `feature_select`. A pattern matches if it occurs anywhere in the name; anchor it with `^`/`$` as needed (e.g. `'^mt-'`).
- **Whole name**: `run_ficture2`, which hands the regexes to `punkst lda4hex` (`std::regex_match`). The pattern must match the entire name, e.g. `'^(mt-.*|Gm[0-9]+)$'` rather than `'^(mt-|Gm[0-9])'`.

---
## Filtering a feature table (`feature_filtering`)

`feature_filtering` filters the rows of a CSV/TSV (e.g., `feature.clean.tsv.gz`) by feature name, using explicit lists, substrings, regex patterns, or a feature-type reference file. The filtered file can then be used downstream, for example as [`run_ficture2 --in-feature-ficture`](./run_ficture2.md#inputoutput-parameters).

### Requirements

- Input CSV/TSV with a header line and a feature-name column (plain text or gzipped).
- A reference file of feature types when using `--include-feature-type-regex`.

### Actions

!!! warning "Action Specifications"
    No action runs by default. Activate at least one filter criterion via the [filter parameters](#filters-parameters).

- Include criteria are restrictive (a feature must satisfy all provided include constraints).
- Exclude criteria are subtractive (a feature matching any exclude constraint is removed).

### Example Usage

#### 1) Filter by Regex

This example demonstrates an exclude-only pattern with `--exclude-feature-regex`; conversely, you could keep only matching features using `--include-feature-regex`.

```bash
cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.filtered.tsv.gz \
  --exclude-feature-regex '^(BLANK|NegCon|NegPrb)'
```

#### 2) Filter by Lists

Include-only and exclude-only lists can be used alone or together. Below is an example using `--include-feature-list`.

```bash
include_list=/path/to/input_keep_genes.txt
cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.filtered.tsv.gz \
  --out-record /path/to/sge/feature.filter_record.tsv \
  --csv-colname-feature-name gene \
  --include-feature-list ${include_list} \
  --log
```

#### 3) Filter by Substrings

Below, include and exclude are shown together for illustration; you may use either one independently.

```bash
cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.neuronal.tsv.gz \
  --include-feature-substr Neur \
  --exclude-feature-substr Pseudo
```

#### 4) Filter by Feature Types

These examples keep features of a given type using a reference file; adapt the regex to the types you want to keep.

```bash
# Example 4A includes `protein_coding` features, selecting columns by name in a TSV reference file.

ref_tsv=/path/to/refs/feature_types.tsv

cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.protein_coding.tsv.gz \
  --include-feature-type-regex '^protein_coding$' \
  --feature-type-ref ${ref_tsv} \
  --feature-type-ref-colname-name gene \
  --feature-type-ref-colname-type type

# Example 4B includes `lncRNA` features, selecting columns by 0-based index in a CSV reference file:

ref_csv=/path/to/refs/gene_types.csv

cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.lncRNA.tsv.gz \
  --include-feature-type-regex '^lncRNA$' \
  --feature-type-ref ${ref_csv} \
  --feature-type-ref-delim "," \
  --feature-type-ref-colidx-name 0 \
  --feature-type-ref-colidx-type 2
```

### Parameters

#### Input/Output Parameters
- `--in-csv` (str, required): Input CSV/TSV with a header line (plain text or gzipped).
- `--out-csv` (str, required): Output CSV/TSV of the rows kept, gzipped (`.gz` is appended if missing).
- `--out-record` (str): Optional table with the `feature` and its `filtering` reason (empty when kept).
- `--csv-colname-feature-name` (str): Feature column name in the input file (default: `gene`).
- `--csv-delim` (str): Delimiter of the input, the output and the record (default: `\t`).
- `--chunksize` (int): Chunk size for streaming reads (default: 50000).
- `--log` (flag): Write logs to `feature_filtering.log` next to the outputs.

#### Filters Parameters

Regex and substring filters match anywhere in the name (substring semantics). List files are plain text with one feature name per line (no header).

##### List Filters Parameters
- `--include-feature-list` (str): Path to a file with feature names to include.
- `--exclude-feature-list` (str): Path to a file with feature names to exclude.

##### Substring Filters Parameters
- `--include-feature-substr` (str): Include features containing the substring.
- `--exclude-feature-substr` (str): Exclude features containing the substring.

##### Regex Filters Parameters
- `--include-feature-regex` (str): Include features matching the regex.
- `--exclude-feature-regex` (str): Exclude features matching the regex.

##### Type Filters Parameters
!!! info "Reference file columns (name vs. index)"
    The reference file is always read with its first line as the header (lines starting with `#` are skipped). Select the feature-name and feature-type columns either by header name (`--feature-type-ref-colname-*`) or by 0-based position (`--feature-type-ref-colidx-*`).

- `--include-feature-type-regex` (str): Include features whose type matches the regex (e.g., `^protein_coding$`).
- `--feature-type-ref` (str): Reference file with feature name and type columns.
- `--feature-type-ref-delim` (str): Delimiter of the reference file (default: `\t`).
- `--feature-type-ref-colname-name` (str): Column name for feature name.
- `--feature-type-ref-colname-type` (str): Column name for feature type.
- `--feature-type-ref-colidx-name` (int): 0-based column index for feature name.
- `--feature-type-ref-colidx-type` (int): 0-based column index for feature type.

### Output

A gzipped CSV/TSV with the rows whose feature passed the filters. If `--out-record` is set, a table mapping each feature to a `filtering` reason; empty when a feature is kept.

---
## Resolving filters into one list (`feature_select`)

`feature_select` resolves include/exclude lists and regexes into a **single feature list file**. It is used:

- by [`run_ficture2_multi`](./run_ficture2_multi.md#feature-filters), which resolves its four filter flags (`--mode include`) against `multi.union_features.tsv` into `multi.selected_features.tsv` and passes it to `lda4hex --features`, restricting the factor analysis without altering the tiled transcripts or the feature list that packaging reads;
- by [`run_ficture2_multi_cells`](./run_ficture2_multi_cells.md#feature-filters), to fold a list and a regex of the same polarity into one list before `spatula`;
- by hand, for [`sge_convert`](./sge_convert.md#feature-filtering-parameters) on the MEX platforms, where a list and a regex of the same polarity cannot be combined.

### Example Usage

Fold an include list and an include regex into one list (the include list is the universe):

```bash
cartloader feature_select \
  --mode include \
  --include-list /path/to/panel_genes.txt \
  --include-regex '^(Cd|Il)' \
  --out /path/to/selected_genes.txt
```

Fold an exclude list and an exclude regex into one list of features to drop:

```bash
cartloader feature_select \
  --mode exclude \
  --in-features /path/to/all_gene_names.txt \
  --exclude-list /path/to/drop_genes.txt \
  --exclude-regex '^(BLANK|NegCon|NegPrb)' \
  --out /path/to/drop_all.txt
```

### Parameters

- `--out` (str, required): Path to the output feature list.
- `--mode` (`include`|`exclude`): `include` (default) writes the features to keep; `exclude` writes the features to drop.
- `--in-features` (str): The universe of features: a TSV whose first column is the feature name and whose optional second column is a count (e.g., `multi.union_features.tsv`); lines starting with `#` are skipped. Required unless `--mode include` with `--include-list`, in which case the include list is the universe.
- `--include-list`, `--exclude-list` (str): Files of feature names to include/exclude (one per line; only the first tab-separated column is read).
- `--include-regex`, `--exclude-regex` (str): Regexes of feature names to include/exclude (substring semantics).
- `--allow-empty` (flag): Write an empty output instead of failing when no feature survives.

### Output

- `--mode include`: the universe features that are in `--include-list` (if given), not in `--exclude-list`, match `--include-regex` (if given), and do not match `--exclude-regex`, in universe order. When every selected feature has a count in the universe, the file starts with a `#feature<TAB>total_count` header and carries the counts; otherwise it lists names only.
- `--mode exclude`: `--exclude-list` united with the `--in-features` entries matching `--exclude-regex`, one name per line. Include flags are not allowed in this mode, and `--exclude-regex` requires `--in-features`.
