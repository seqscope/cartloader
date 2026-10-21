# Step-by-step

[`run_together`](../reference/run_together.md) chains the `CartLoader` commands for you. This section is for running those commands yourself, one stage at a time.

!!! tip "When to run the individual commands"
    - To learn what each stage does and what it produces.
    - To customize a stage beyond what `run_together` exposes (for example a parameter of `run_ficture2` that has no `run_together` option).
    - To rerun or replace a single stage, or to combine `CartLoader` with your own tools.

    For a standard run, use `run_together`: it applies platform defaults, runs the stages in the right order as one resumable Makefile, and packages the result for CartoScope.

## The Stages and Their Commands

| Stage | Command(s) | What it does |
|---|---|---|
| SGE conversion | [`sge_convert`](../reference/sge_convert.md), [`reformat_cosmx`](../reference/reformat_cosmx.md) | Converts raw platform output into the unified transcript format. |
| SGE stitching and orientation | [`sge_stitch`](../reference/sge_stitch.md), [`sge_orientate`](../reference/sge_orientate.md) | Combines tiles or sections into one SGE; rotates or flips coordinates. |
| Spatial factors (FICTURE) | [`run_ficture2`](../reference/run_ficture2.md), [`run_ficture2_multi`](../reference/run_ficture2_multi.md) | Trains or projects `FICTURE` models and decodes pixels; `_multi` trains one model across samples. |
| Cell analyses | [`run_ficture2_multi_cells`](../reference/run_ficture2_multi_cells.md), [`import_cell`](../reference/import_cell.md), [`import_square`](../reference/import_square.md) | Analyzes segmented cells, or imports platform cell and bin analyses (Xenium Ranger, Space Ranger). |
| Background images | [`import_image`](../reference/import_image.md) | Converts DAPI, H&E, and other images into aligned PMTiles. |
| Packaging | [`run_cartload2`](../reference/run_cartload2.md), [`run_cartload2_multi`](../reference/run_cartload2_multi.md) | Packages the SGE, factors, and cell layers into PMTiles with a `catalog.yaml`. |
| Publishing | [`upload_aws`](../reference/upload_aws.md), [`upload_zenodo`](../reference/upload_zenodo.md) | Uploads a packaged directory to S3 or Zenodo. |
| AI annotation | [`anno_cartload_folder`](../reference/anno_cartload_folder.md), [`annotate_factors_with_llm`](../reference/annotate_factors_with_llm.md) | Labels factors and cell clusters with an LLM. |

`run_together` uses the multi-sample commands (`run_ficture2_multi`, `run_ficture2_multi_cells`) even for a single sample, so its outputs follow their layout. The single-sample `run_ficture2` and `run_cartload2` remain available for step-by-step use.

## Walkthroughs

- **Per platform**: each [platform tutorial](./intro.md#platform-tutorials) has a "step by step" part that runs `sge_convert` → `run_ficture2` → `run_cartload2` → upload on a small dataset, after the equivalent `run_together` command.
- **Multiple samples**: the [Human Cortex Multi-sample Tutorial](./multi_sample.md) runs `run_ficture2_multi` and `run_cartload2_multi` on five human cortex MERFISH samples.

!!! question
    - [Which interface should I use?](../faq/choose_interface.md)
    - [How do I resume or restart a run?](../faq/restart_or_resume.md)
