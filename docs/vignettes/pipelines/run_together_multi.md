# End-to-End with `run_together` (multi-sample)

When you have several sections that should share **one FICTURE model**, list them in a sample sheet and `run_together` will:

1. ingest every sample,
2. train a **single joint model** on all of them together,
3. package each sample independently against that shared model.

The dependency graph fans in on the joint model, then fans back out per sample — so `make -j` ingests all samples in parallel, blocks once on the shared model, then packages each sample in parallel:

```text
ingest (sample 1) ─┐
                   ├─► ficture (joint model) ─┬─► cartload 1 ─► images 1 ─► publish 1
ingest (sample 2) ─┘                          └─► cartload 2 ─► images 2 ─► publish 2
```

!!! info "When to use a joint model"
    Sharing one model across sections makes factors directly comparable between samples and is the recommended setup for replicates or a cohort. For unrelated samples, run them independently (see [batches](#batches-of-independent-samples)). Background: [Why multi-sample?](../../faq/why_multi-sample.md)

---
## 1. List the samples in a TSV sheet

A multi-sample run is just a single-sample run with `--samples samples.tsv` in place of `--in-dir`. The sheet is a wide table of input roles; for the common case you only need `id` + `in_dir`:

```
id      in_dir
rep1    /data/kidney/rep1/outs
rep2    /data/kidney/rep2/outs
rep3    /data/kidney/rep3/outs
```

## 2. Run

```bash
cartloader run_together --platform 10x_xenium \
    --samples samples.tsv --out-dir /path/to/output/kidney-cohort \
    --width 18 --n-factor 24 -j 10
```

All rows share one `--out-dir`, so they train **one joint model**. Everything else (FICTURE mode, images, cell analyses) is assembled exactly as in the single-sample run and applied to every sample. Add `--dry-run` to inspect `run_together.mk` and `run_together.resolved.json` first.

!!! info "Output layout"
    Packaging is delegated to [`run_cartload2_multi`](../../reference/run_cartload2_multi.md): each sample is written to a self-contained `cartl/<multi_id>-<sample_id>/` directory (with `<multi_id>` defaulting to the `--out-dir` basename), plus a `cartl/multi-catalog.yaml` that links every per-sample `catalog.yaml` and copies the shared factor files (`post`/`rgb`/`de`/`info`/`umap`) into the `cartl/` root as a unified `factors:` map. Upload the whole `cartl/` directory to S3 as one unit.

!!! tip "Extra role columns"
    The sheet can also carry per-sample inputs — `raw_transcript` (a raw CSV to ingest), `transcript` (pre-converted TSV, skips ingest), `xy`, `boundaries`, `clusters`, `mex`, `cellxgene` (a cell×gene matrix CSV, e.g. MERSCOPE `cell_by_gene.csv`):
    ```
    id    in_dir       transcript                     boundaries               cellxgene
    rep1  /data/rep1
    rep2               /data/rep2/transcripts.tsv.gz  /data/rep2/bounds.csv.gz
    rep3  /data/rep3                                                            /data/rep3/cell_by_gene.csv
    ```
    See the [reference → sample sheet](../../reference/run_together_inputs.md#the-sample-sheet-is-a-wide-table-of-input-roles) for the full column list, aliases, and the three transcript sources.

### When to use JSON instead

Escalate to a `--config` JSON when samples need **different** settings, or to add analyses/images beyond the profile defaults. JSON augments the CLI (append-by-default); a `samples` block gives full per-sample control:

```jsonc
{
  "samples": [
    { "id": "rep3", "in_dir": "/data/kidney/rep3/outs",
      "hne": "/data/kidney/rep3/he.tif",
      "cartload": { "bin_count": 300 } }
  ]
}
```

See the [configuration reference](../../reference/run_together_inputs.md#mode-3-full-config-json) for every key.

---
## Visium HD example

The same shape works for Visium HD; the `10x_visium_hd` profile adds the square/cell imports and (with a per-sample `hne` column) the H&E layer automatically. Add `hne` as a sheet column:

```
id      in_dir             hne
38088   /data/SI_38088/outs  /data/HE_38088.tif
39685   /data/SI_39685/outs  /data/HE_39685.tif
```

```bash
cartloader run_together --platform 10x_visium_hd \
    --samples samples.tsv --out-dir /path/to/output/prostate-cohort -j 10
```

---
## Publishing (opt-in)

Publishing is CLI-driven — no config block. Enable annotation (`--anno`, needs `--tissue`/`--organism`) and/or S3 upload (`--s3-upload`) for every sample:

```bash
cartloader run_together --platform 10x_xenium \
    --samples samples.tsv --out-dir OUT --width 18 --n-factor 24 -j 10 \
    --anno --tissue "Kidney" --organism human \
    --s3-upload --collection kidney-cohort
```

Each sample uploads to `s3://cartostore/data/batch=<YYYY_MM>/<collection>/<multi_id>-<sample_id>/`, where `<YYYY_MM>` defaults to the current date (override `--batch`) and `<collection>` defaults to the out-dir basename (override `--collection`). See the [reference → Publish](../../reference/run_together.md#publishing) for all flags (`--s3-prefix`, `--aws-profile`, annotation model/threads, …).

---
## Batches of independent samples

For many **unrelated** samples (each its own model), use `--out-root` instead of `--out-dir`:

```bash
cartloader run_together --platform 10x_xenium \
    --samples samples.tsv --out-root /path/to/output/2026_07 -j 6
```

Each row becomes its own `<out_root>/<id>/` with an independent model. See [Specifying Inputs → Multi-sample](../../reference/run_together_inputs.md#mode-2-multi-sample-sample-sheet).

---
## Next steps

- Full option and configuration reference: [`run_together` reference](../../reference/run_together.md).
- Supported platforms and expected inputs: [Supported Platforms & Inputs](../../reference/supported_platforms.md).
- A hand-built multi-sample walkthrough (lower-level modules): [Human Cortex Multi-sample Tutorial](../multi_sample.md).
