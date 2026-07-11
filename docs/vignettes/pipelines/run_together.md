# End-to-End with `run_together` (single sample)

`run_together` runs a **complete CartoScope pipeline with a single command** — ingest → FICTURE → cell decode → asset packaging → image import → (optional) publish. It builds one Makefile and runs it, so the pipeline is **resumable** and **parallel** out of the box.

This tutorial mirrors the [step-by-step Xenium tutorial](./xenium.md) but replaces the whole manual chain with one command, using the same public 10x Xenium human lung cancer dataset.

!!! info "Which entry point should I use?"
    `run_together` is one of several ways to drive CartLoader. It shines for **multi-platform** and **multi-sample** runs and for one-command convenience. The per-platform orchestrators [`run_xenium`](../../reference/run_xenium.md) / [`run_visiumhd`](../../reference/run_visiumhd.md) and the individual modules remain fully supported. See [Which interface should I use?](../../faq/choose_interface.md).

---
## Prepare input

Download and unpack the dataset exactly as in the [step-by-step Xenium tutorial → Prepare Input](./xenium.md#prepare-input). After unzipping, `${work_dir}/raw` is the Xenium Ranger output directory.

The `10x_xenium` **profile** reads this standard layout automatically — you do **not** list these files yourself:

```
${work_dir}/raw/
├── transcripts.parquet                 # or transcripts.csv.gz
├── cell_boundaries.csv.gz              # cell decode
├── cells.csv.gz                        # cell centroids (x_centroid, y_centroid)
├── analysis/clustering/gene_expression_graphclust/clusters.csv   # → "xeniumranger" prefix
└── morphology_focus/morphology_focus_000{0,1,2,3}.ome.tif        # dapi/boundary/rna/protein
```

See the full [expected-input reference](../../reference/supported_platforms.md#10x_xenium) for every path the profile looks for.

---
## Define ID and parameters

```bash
DATA_ID="xenium-v1-humanlung-cancer-ffpe"   # dataset name; also names the output folder
train_width=18                              # FICTURE hexagon width (µm)
n_factor=24                                 # number of factors (comma-separated for several)
n_jobs=10
```

!!! warning "External tools on PATH"
    `run_together` delegates to CartLoader modules, which expect `spatula`, `punkst` (FICTURE2), `tippecanoe`, `go-pmtiles`, `gdal`, and `pigz` available (on `PATH`, or as built repo submodules for `spatula`/`punkst`). Unlike `run_xenium`, `run_together` does not take per-tool path flags.

---
## Run the pipeline

=== "Locally"

    **Set up the environment**
    {%
    include-markdown "../../../includes/includemd_vigenettes_setupenv.md"
    %}

    **Command**
    ```bash
    cartloader run_together \
        --platform 10x_xenium \
        --in-dir  ${work_dir}/raw \
        --out-dir ${work_dir}/output/${DATA_ID} \
        --width ${train_width} \
        --n-factor ${n_factor} \
        -j ${n_jobs} --threads ${n_jobs}
    ```

=== "Via Docker"

    **Set up the environment**
    {%
    include-markdown "../../../includes/includemd_vigenettes_setupenv_docker.md"
    %}

    **Command**
    ```bash
    docker run -it --rm \
        -v ${work_dir}:/data \
        weiqiuc/cartloader:${docker_tag} \
        run_together \
        --platform 10x_xenium \
        --in-dir  /data/raw \
        --out-dir /data/output/${DATA_ID} \
        --width ${train_width} \
        --n-factor ${n_factor} \
        -j ${n_jobs} --threads ${n_jobs}
    ```

!!! tip "Preview before running"
    Add `--dry-run` to write the Makefile and print every command (`make -n`) without executing. Inspect `${work_dir}/output/${DATA_ID}/run_together.mk` and `run_together.resolved.json` to see exactly what will run.

### Projection-only mode

To reuse models from a **previous FICTURE run** instead of training new ones, point `--project-models` at that run's FICTURE directory. `run_together` reads its `ficture.params.json` and re-projects every model it lists — no LDA training runs.

```bash
cartloader run_together --platform 10x_xenium \
    --in-dir ${work_dir}/raw --out-dir ${work_dir}/output/${DATA_ID} \
    --project-models /path/to/previous/fic --width 12 -j ${n_jobs}
```

---
## Output

```
${work_dir}/output/${DATA_ID}/
├── run_together.mk               # the generated pipeline (Makefile)
├── run_together.resolved.json    # fully-resolved settings (provenance)
├── mk/                           # per-stage flag files that drive make
├── tsv/<id>/transcripts.unsorted.tsv.gz
├── fic/                          # FICTURE results (per sample under fic/samples/<id>/)
└── cartl/<id>/                   # packaged PMTiles + catalog.yaml  ← deploy this
```

The `cartl/<id>/catalog.yaml` plus its PMTiles is the deployable CartoScope asset. (A joint multi-sample run instead produces `cartl/<multi_id>-<sample_id>/` per sample plus a `cartl/multi-catalog.yaml` — see the [multi-sample tutorial](./run_together_multi.md).) See the per-module output details in [`sge_convert`](../../reference/sge_convert.md#output), [`run_ficture2`](../../reference/run_ficture2.md#output), and [`run_cartload2`](../../reference/run_cartload2.md#output).

---
## Resume, re-run, and publish

- **Resume after a failure:** just re-run the same command (or `make -f .../run_together.mk -j N`). Completed stages are skipped via their flag files.
- **Run part of the pipeline:** `--only ingest,ficture` or `--skip images`. Excluded upstream stages are assumed already done.
- **Force a clean rebuild:** `--restart` (`make -B`).
- **Publish (opt-in):** annotation + S3 upload run only with a `publish` block in a `--config` file **and** the `--publish` flag. See the [multi-sample tutorial](./run_together_multi.md) and the [reference](../../reference/run_together.md#publishing-opt-in).

---
## Next steps

- Process several sections under **one joint FICTURE model**: [End-to-End with `run_together` (multi-sample)](./run_together_multi.md).
- Full option and configuration reference: [`run_together` reference](../../reference/run_together.md).
- Supported platforms and the exact inputs each expects: [Supported Platforms & Inputs](../../reference/supported_platforms.md).
