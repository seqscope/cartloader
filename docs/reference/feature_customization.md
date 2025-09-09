# Feature Filtering (Add-on)

## Overview

While `run_ficture2` allows filtering by regex, CartLoader also provides `feature_filtering`, which filters features (e.g., genes) from a CSV/TSV using explicit lists, substrings, regex patterns, or a type reference file. It writes a filtered CSV/TSV (gzipped) and can optionally export an annotation record describing which features were removed and why. The curated feature file can then be used downstream, including FICTURE analysis.

## Example Usage

### 1) Filter by Regex

This example demonstrates an exclude-only pattern with `--exclude-feature-regex`; conversely, you could keep only matching features using `--include-feature-regex`.

```bash
cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.filtered.tsv.gz \
  --exclude-feature-regex '^(BLANK.*$|NegCon.*$|NegPrb.*$)'
```

### 2) Filter by Lists

CartLoader supports include-only and exclude-only list-based filtering; pick the one that fits your use case (or combine them). Below is an example using `--include-feature-list`.

```bash
include_list=/path/to/input_keep_genes.txt 
cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.filtered.tsv.gz \
  --out-record /path/to/sge/feature.filter_record.tsv \
  --csv-colname-feature-name gene \
  --include-feature-list ${include_list}\
  --log
```

### 3) Filter by Substrings

Below, include and exclude are shown together for illustration; you may use either one independently.

```bash
cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.neuronal.tsv.gz \
  --include-feature-substr Neur \
  --exclude-feature-substr Pseudo
```

### 4) Filter by Feature Types

These examples show include-by-type filtering using a reference file; adapt the regex to the types you want to keep. Below provides two examples:

```bash
# Example 4A includes `protein_coding` features and uses column names and a TSV reference file.

ref_tsv=/path/to/refs/feature_types.tsv

cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.protein_coding.tsv.gz \
  --include-feature-type-regex '^protein_coding$' \
  --feature-type-ref ${ref_tsv} \
  --feature-type-ref-colname-name gene \
  --feature-type-ref-colname-type type

# Example 4B includes `lncRNA` features using column indices (0-based) and a CSV reference file:

ref_csv=/path/to/refs/gene_types.csv 

cartloader feature_filtering \
  --in-csv /path/to/sge/feature.clean.tsv.gz \
  --out-csv /path/to/sge/feature.lncRNA.tsv.gz \
  --include-feature-type-regex '^lncRNA$' \
  --feature-type-ref ${ref_csv} \
  --csv-delim "," \
  --feature-type-ref-colidx-name 0 \
  --feature-type-ref-colidx-type 2
```

## Parameters

!!! warning "Filtering Action"
    At least one filtering criterion is required; otherwise the command errors.

    Also, please note that:

    - Include criteria are restrictive (a feature must satisfy all provided include constraints).
    - Exclude criteria are subtractive (a feature matching any exclude constraint is removed).

### Input/Output
- `--in-csv` (str, required): Input CSV/TSV of features (plain text or gzipped).
- `--out-csv` (str, required): Output filtered CSV/TSV; gzipped. If extension is not `.gz`, the tool appends `.gz`.
- `--out-record` (str): Optional TSV with `feature` and `filtering` reason per feature.
- `--csv-colname-feature-name` (str): Feature column name in the input file (default: `gene`).
- `--csv-delim` (str): Field delimiter for the input (applies to output as well; default: `\t`).
- `--chunksize` (int): Chunk size for streaming reads (default: 50000).
- `--log` (flag): Write logs to `feature_filtering.log` next to outputs.

### List Filters
- `--include-feature-list` (str): Path to a file with feature names to include.
- `--exclude-feature-list` (str): Path to a file with feature names to exclude.

### Substring Filters
- `--include-feature-substr` (str): Include features containing the substring.
- `--exclude-feature-substr` (str): Exclude features containing the substring.

### Regex Filters
- `--include-feature-regex` (str): Include features matching the regex.
- `--exclude-feature-regex` (str): Exclude features matching the regex.

### Type Filters
If the reference file has a header row, specify the feature-name and feature-type columns with `--feature-type-ref-colname-*`; otherwise, use the 0-based index flags `--feature-type-ref-colidx-*`.

- `--include-feature-type-regex` (str): Include by feature type (e.g., `^protein_coding$`).
- `--feature-type-ref` (str): Reference file with feature name and type columns.
- `--feature-type-ref-delim` (str): Delimiter for the reference file (default: `\t`).
- `--feature-type-ref-colname-name` (str): Column name for feature name.
- `--feature-type-ref-colname-type` (str): Column name for feature type.
- `--feature-type-ref-colidx-name` (int): 0-based column index for feature name.
- `--feature-type-ref-colidx-type` (int): 0-based column index for feature type.

## Outputs
- `--out-csv`: Filtered rows from the input CSV/TSV, gzipped. Overwrites existing file.
- `--out-record` (optional): TSV mapping each feature to a `filtering` reason; empty when a feature is kept.
