# Generic Platform Tutorial

Use the [generic profile](../../reference/platforms/generic.md) (`--platform generic`) for any platform without a dedicated profile, or whenever your data is a plain transcript table: one row per molecule (or per spatial bin and gene) with X, Y, gene, and optionally a count and a cell ID.

To show the workflow on real data, this tutorial treats the [Xenium starter](./xenium.md) table, a small region of the adult mouse brain, as if it came from an unknown platform. It has the typical quirks of a table from another tool: comma-separated, its own column names, no count column, and a cell-ID column. (For actual Xenium data, use the [Xenium profile](../../reference/platforms/xenium.md) instead.)

It first runs the whole analysis with a single [`run_together`](../../reference/run_together.md) command, then walks through the same analysis step by step.

---

## Run with `run_together`

### Input Data

Download the example data, hosted on Zenodo ([10.5281/zenodo.17953582](https://doi.org/10.5281/zenodo.17953582)), into a work directory:

```bash
work_dir=/path/to/work/directory        # replace with your work directory
mkdir -p ${work_dir} && cd ${work_dir}
wget https://zenodo.org/records/17953582/files/xenium_starter.raw.tar.gz
tar --strip-components=1 -zxvf xenium_starter.raw.tar.gz   # extracts input.tsv.gz
```

`input.tsv.gz` is comma-separated despite its name. Its first lines:

```
"transcript_id","cell_id","overlaps_nucleus","feature_name","x_location","y_location","z_location","qv"
281827164036151,133793,1,"Sox10",2350.0232,4153.6846,16.592316,40.0
```

### Describe the Table

Tell `run_together` how to read the table:

| What `run_together` needs | In this file | Option |
|---|---|---|
| X and Y coordinates | `x_location`, `y_location` | `--colname-transcript-x`, `--colname-transcript-y` |
| Gene name | `feature_name` | `--colname-transcript-feature` |
| Count | none: each row is one molecule | `--colname-transcript-count` (omitted here) |
| Coordinate unit | micrometers | `--units-per-um 1` (input units per µm) |
| Cell ID (optional) | `cell_id` | `--colname-transcript-cell cell_id` and `--colname-cell cell_id` |
| Delimiter | comma | `ingest.extra_flags` in a config file (below); tab is the default |

The delimiter has no command-line option, so set it in a small config file:

```bash
cat > generic_starter.json <<'EOF'
{ "ingest": { "extra_flags": ["--csv-delim ,"] } }
EOF
```

### Run

Run the pipeline, either with a local [installation](../../installation.md) or with the [Docker image](../quickstart/run_in_docker.md):

=== "Locally"

    ```bash
    cartloader run_together \
      --platform generic \
      --raw-transcript ./input.tsv.gz \
      --config generic_starter.json \
      --colname-transcript-x x_location \
      --colname-transcript-y y_location \
      --colname-transcript-feature feature_name \
      --colname-transcript-cell cell_id \
      --colname-cell cell_id \
      --units-per-um 1 \
      --out-dir ./out/generic_hippo \
      --id generic_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

=== "Docker"

    ```bash
    docker run --rm -v ${work_dir}:/data hyunminkang/cartloader:latest \
      run_together \
      --platform generic \
      --raw-transcript /data/input.tsv.gz \
      --config /data/generic_starter.json \
      --colname-transcript-x x_location \
      --colname-transcript-y y_location \
      --colname-transcript-feature feature_name \
      --colname-transcript-cell cell_id \
      --colname-cell cell_id \
      --units-per-um 1 \
      --out-dir /data/out/generic_hippo \
      --id generic_hippo \
      --width 12 \
      --n-factor 6,12 \
      --threads 4 -j 4
    ```

What happens:

- **SGE conversion** reads the table with the column mapping above and writes the unified format. The generic profile drops negative-control and unassigned probes by default (`^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated|System|NCS-|NCP-)`); set `--ingest-exclude-feature-regex` to change that.
- **FICTURE** trains one model per `--width` × `--n-factor` combination and decodes every pixel.
- **Cell analysis** runs because `--colname-cell` is given: transcripts are grouped by `cell_id` into cells, which are clustered with the largest model. Transcripts outside any cell (`cell_id` of `-1` here) are ignored. Without the two cell options, the run stays pixel-level.
- **Packaging** writes the deployable result to `out/generic_hippo/cartl/generic_hippo/`: a `catalog.yaml` and the PMTiles it lists.

Add `--dry-run` to print the commands without running them.

### Adapting It to Your Data

- **A tab-separated table** needs no config file: tab is the default delimiter.
- **A count column** (e.g. UMIs per spatial bin): name it with `--colname-transcript-count`.
- **Other units**: set `--units-per-um` to the number of coordinate units per micrometer, e.g. `1000` for nanometers (see [How do I define scaling?](../../faq/define_units_per_um.md)).
- **No cell IDs**: omit `--colname-transcript-cell` and `--colname-cell`.
- **Cell centroids or polygons** from your segmentation can be added with `--in-cell-xy` and `--in-cell-boundary`; see the [generic platform page](../../reference/platforms/generic.md) for their formats, for several samples at once, and for a full config.

{%
  include-markdown "../../../includes/includemd_vigenettes_run_together_diff.md"
%}

---

## Step by Step (Individual Commands)

The same pixel-level analysis, one command at a time. (The cell analysis that `run_together` adds is run by [`run_ficture2_multi_cells`](../../reference/run_ficture2_multi_cells.md); it is easiest to get through `run_together`.)

### Set Up the Environment

{%
  include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
%}

### Input Data

Use the data downloaded [above](#input-data).

### Define ID and Parameters

```bash
# Unique identifier for your dataset
DATA_ID="generic_hippo"                  # change this to reflect your dataset name

# LDA parameters
train_width=12                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

### SGE Format Conversion

Convert the table to the unified SGE format with [`sge_convert`](../../reference/sge_convert.md) in generic mode. The `--csv-*` options play the role of `run_together`'s `--colname-transcript-*` options:

```bash
cartloader sge_convert \
  --makefn sge_convert.mk \
  --platform generic \
  --in-csv ./input.tsv.gz \
  --csv-delim , \
  --csv-colname-x x_location \
  --csv-colname-y y_location \
  --csv-colname-feature-name feature_name \
  --units-per-um 1 \
  --out-dir ./sge \
  --exclude-feature-regex '^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated|System|NCS-|NCP-)' \
  --sge-visual \
  --n-jobs ${n_jobs}
```

Without `--csv-colname-count`, each row counts as one molecule.

### `FICTURE` Analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
%}

### `CartLoader` Asset Packaging

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
%}

### Upload to Data Repository

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md" preserve-includer-indent=false
%}

---

## Output Data

The `run_together` run writes everything under `out/generic_hippo/` (see [Where are outputs written?](../../faq/find_outputs.md)); the deployable result is `out/generic_hippo/cartl/generic_hippo/`. For the outputs of each individual command, see the reference pages for [`sge_convert`](../../reference/sge_convert.md), [`run_ficture2`](../../reference/run_ficture2.md), and [`run_cartload2`](../../reference/run_cartload2.md).
