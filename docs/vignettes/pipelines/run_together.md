# Single Sample with `run_together`

`run_together` runs a **complete CartoScope pipeline with a single command**: SGE conversion → FICTURE → cell analyses → background images → asset packaging → (optional) publishing. It builds one Makefile and runs it, so the pipeline is **resumable** and **parallel** out of the box.

This guide processes a public 10x Xenium human lung cancer dataset (a full section, larger than the [quick start](../quickstart/intro.md) examples).

!!! info "The recommended entry point"
    `run_together` is the recommended way to run `CartLoader` on any supported platform, for one sample or [many](./run_together_multi.md). The individual commands remain available for [step-by-step](../step_by_step.md) use; the older `run_xenium` / `run_visiumhd` orchestrators are [legacy](../../reference/legacy.md). See [Which interface should I use?](../../faq/choose_interface.md).

---
## Prepare Input

Download the Xenium Ranger output from the [10x Genomics Dataset portal](https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard) and unzip it:

```bash
work_dir=/path/to/work/directory
mkdir -p ${work_dir}/raw && cd ${work_dir}/raw

wget https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_humanLung_Cancer_FFPE/Xenium_V1_humanLung_Cancer_FFPE_outs.zip
# or: curl -O https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_humanLung_Cancer_FFPE/Xenium_V1_humanLung_Cancer_FFPE_outs.zip
unzip Xenium_V1_humanLung_Cancer_FFPE_outs.zip
```

`${work_dir}/raw` is now the Xenium Ranger output directory. The `10x_xenium` **profile** reads its standard layout automatically, so you do **not** list these files yourself:

```
${work_dir}/raw/
├── transcripts.parquet                 # or transcripts.csv.gz
├── cell_boundaries.csv.gz              # cell boundaries for the cell analyses
├── cells.csv.gz                        # cell centroids (x_centroid, y_centroid)
├── analysis/clustering/gene_expression_graphclust/clusters.csv   # Xenium Ranger clusters
└── morphology_focus/morphology_focus_000{0,1,2,3}.ome.tif        # dapi/boundary/rna/protein images
```

See the [Xenium platform page](../../reference/platforms/xenium.md) for every path the profile looks for.

---
## Define ID and Parameters

```bash
DATA_ID="xenium-v1-humanlung-cancer-ffpe"   # dataset name; names the output folder and the packaged assets
train_width=18                              # FICTURE hexagon width (µm)
n_factor=24                                 # number of factors (comma-separated for several)
n_jobs=10
```

---
## Run the Pipeline

=== "Locally"

    Requires a local [installation](../../installation.md), including the tools built by `submodules/build.sh`; `run_together` finds them there, so it takes no tool-path options.

    ```bash
    cartloader run_together \
        --platform 10x_xenium \
        --in-dir  ${work_dir}/raw \
        --out-dir ${work_dir}/output/${DATA_ID} \
        --id ${DATA_ID} \
        --width ${train_width} \
        --n-factor ${n_factor} \
        -j ${n_jobs} --threads ${n_jobs}
    ```

=== "Via Docker"

    Uses the prebuilt image, which includes every tool (see [Quick Start with Docker](../quickstart/run_in_docker.md)).

    ```bash
    docker run --rm \
        -v ${work_dir}:/data \
        hyunminkang/cartloader:latest \
        run_together \
        --platform 10x_xenium \
        --in-dir  /data/raw \
        --out-dir /data/output/${DATA_ID} \
        --id ${DATA_ID} \
        --width ${train_width} \
        --n-factor ${n_factor} \
        -j ${n_jobs} --threads ${n_jobs}
    ```

With this profile, one run produces: the converted SGE; one `FICTURE` model per `--width` × `--n-factor` combination with pixel-level decoding; two cell analyses (`cartloader`, a `FICTURE` analysis of the segmented cells, and `xeniumranger`, an import of the Xenium Ranger clusters); the four morphology images; and the packaged PMTiles.

!!! tip "Preview before running"
    Add `--dry-run` to write the Makefile and print every command (`make -n`) without executing. Inspect `${work_dir}/output/${DATA_ID}/run_together.mk` and `run_together.resolved.json` to see exactly what will run.

### Projection-only Mode

To reuse models from a **previous FICTURE run** instead of training new ones, point `--project-models` at that run's FICTURE directory. `run_together` reads its `ficture.params.json` and re-projects every model it lists; no LDA training runs.

```bash
cartloader run_together --platform 10x_xenium \
    --in-dir ${work_dir}/raw --out-dir ${work_dir}/output/${DATA_ID} --id ${DATA_ID} \
    --project-models /path/to/previous/fic --width 12 -j ${n_jobs}
```

---
## Output

```
${work_dir}/output/${DATA_ID}/
├── run_together.mk               # the generated pipeline (Makefile)
├── run_together.resolved.json    # fully resolved settings (provenance)
├── mk/                           # per-stage flag files that drive make
├── tsv/<id>/                     # converted SGE (transcripts.unsorted.tsv.gz, features, coordinate ranges)
├── fic/                          # FICTURE results (per sample under fic/samples/<id>/)
└── cartl/<id>/                   # packaged PMTiles + catalog.yaml  ← deploy this
```

The packaged directory is named `<out-dir name>-<id>`, shortened to `<id>` when the two are equal, as in the command above. Without `--id`, the sample ID is `rep1`, so the directory would be `cartl/xenium-v1-humanlung-cancer-ffpe-rep1/`.

The `catalog.yaml` plus its PMTiles is the deployable CartoScope asset. A joint multi-sample run produces one such directory per sample plus a `cartl/multi-catalog.yaml`; see [Multiple Samples with run_together](./run_together_multi.md). For per-stage output details, see [`sge_convert`](../../reference/sge_convert.md), [`run_ficture2_multi`](../../reference/run_ficture2_multi.md), and [`run_cartload2`](../../reference/run_cartload2.md).

---
## Resume, Re-run, and Publish

- **Resume after a failure:** re-run the same command (or `make -f .../run_together.mk -j N`). Completed stages are skipped via their flag files.
- **Run part of the pipeline:** `--only ingest,ficture` or `--skip images`. Excluded upstream stages are assumed already done.
- **Force a clean rebuild:** `--restart` (`make -B`).
- **Publish (opt-in):** AI annotation and S3 upload are separate flags: `--anno` or `--anno-deep` (with `--tissue`/`--organism`) and `--s3-upload`. See [Publishing](../../reference/run_together.md#publishing).

---
## Next Steps

- Process several sections under **one joint FICTURE model**: [Multiple Samples with run_together](./run_together_multi.md).
- Every option and configuration key: [`run_together` reference](../../reference/run_together.md).
- The exact inputs each platform expects: [Platforms](../../reference/supported_platforms.md).
