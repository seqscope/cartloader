# Platform: 10x Xenium (`--platform 10x_xenium`)

Built-in profile. Ingest runs `sge_convert`. Point `--in-dir` at a **Xenium Ranger output directory** and the profile auto-detects everything below.

## Expected inputs (under `--in-dir`)

First match wins where several paths are listed; optional files are skipped when absent.

| Purpose | Path under `--in-dir` |
|---------|-----------------------|
| Transcripts | `transcripts.csv.gz`, or `transcripts.parquet`, or `transcripts/transcripts.parquet` |
| Cell boundaries | `cell_boundaries.csv.gz` |
| Cell centroids | `cells.csv.gz` (columns `x_centroid`, `y_centroid`) |
| Cluster labels | `analysis/clustering/gene_expression_graphclust/clusters.csv` → the `xeniumranger` factor |
| Morphology images | `morphology_focus/morphology_focus_000{0,1,2,3}.ome.tif` → `dapi`/`boundary`/`rna`/`protein`; or a single `morphology_focus.ome.tif` / `morphology.ome.tif` → `dapi` |

**Defaults:** FICTURE `width=12`, `n_factor=12,24,48`, single-molecule mode; packaging with `--use-pmpoint --bin-count 500`. Two cell analyses are wired: `cartloader` (jointly-decoded factors from boundaries + xy) and `xeniumranger` (the Ranger `clusters.csv`).

---
## 1. Single sample (direct CLI)

```bash
cartloader run_together --platform 10x_xenium \
    --in-dir /data/xenium_lung --out-dir OUT \
    --width 18 --n-factor 24 -j 8 --threads 8
```

Morphology images are auto-detected — nothing to list. Everything (transcripts, boundaries, centroids, clusters, images) comes from the directory.

---
## 2. Multi-sample (sample sheet)

One joint FICTURE model across sections sharing `--out-dir`. For standard Ranger directories you only need `id` + `in_dir`:

```bash
cartloader run_together --platform 10x_xenium \
    --samples samples.tsv --out-dir OUT --width 18 --n-factor 24 -j 10
```

```
id      in_dir
rep1    /data/xenium/rep1/outs
rep2    /data/xenium/rep2/outs
```

**GEO-style layouts** (non-standard filenames / scattered paths) — name each file per role with sheet columns. The raw transcript goes in `raw_transcript` (it still needs ingest; the `transcript` column is for an *already-ingested* TSV and skips ingest). Role paths are forwarded to the per-sample cell import as `--csv-*` overrides:

```
id      raw_transcript                xy                     boundaries               clusters                dapi
rep1    /geo/rep1_transcripts.csv.gz  /geo/rep1_cells.csv.gz /geo/rep1_bounds.csv.gz  /geo/rep1_clusters.csv  /geo/rep1_morphology.ome.tif
rep2    /geo/rep2_transcripts.parquet /geo/rep2_cells.parquet /geo/rep2_bounds.parquet -                       /geo/rep2_morphology_focus.ome.tif
```

- **Parquet is fine.** `transcripts.parquet`, `cells.parquet` and `cell_boundaries.parquet` are converted to `.csv.gz` (`OUT/tsv/<id>/parquet2csv/`) as the first ingest step, with a `NOTE:` at planning time; the run aborts if a conversion fails.
- **No clusters?** The `xeniumranger` import (Ranger's own clusters) needs `clusters`; a sample without it — as a column, or as `analysis/clustering/gene_expression_graphclust/clusters.csv` under `in_dir` — simply skips that import (a `NOTE:` says so). The `cartloader` analysis (clustering recomputed on the shared SGE) needs only `xy` + `boundaries` and still runs.
- **DAPI.** A per-sample `dapi` column takes the image directly. a file ending in `morphology.ome.tif` (any prefix) is the multi-page z-stack and is imported with `--use-middle-page` automatically; `morphology_focus.ome.tif` is the 2D projection and needs nothing special. See [Image Modalities](../run_together_images.md).

---
## 3. Full config (JSON)

Use `--config` when samples differ or to add analyses/images. Example: de-novo base on the CLI, plus a projection model and an extra protein image.

```bash
cartloader run_together --platform 10x_xenium --samples samples.tsv --out-dir OUT \
    --width 18 --n-factor 24 --config extra.json -j 10
```
```jsonc
// extra.json
{
  "ficture": [ { "id": "ref", "mode": "project", "model": "/models/ref.tsv", "width": 12 } ],
  "images":  [ { "id": "cd3", "source": "cd3.ome.tif", "kind": "single", "color": "FF0000" } ]
}
```

---
## Notes

- **Xenium Ranger clusters on joint runs.** The `xeniumranger` analysis relies on **sample-specific** cluster labels, so it cannot be jointly decoded: a single-sample run decodes it via `run_ficture2_multi_cells`, but a **joint run** imports each sample's clusters per-sample via `import_xenium_cell` (appended to that sample's catalog). The `cartloader` factor, whose clustering is recomputed on the shared SGE, is decoded jointly in both cases.
- **Images:** the `morphology_focus_000{0..3}` channels map to `dapi`/`boundary`/`rna`/`protein`. To recolor or add channels, see [Image Modalities](../run_together_images.md).

## See also

- [Specifying Inputs](../run_together_inputs.md) · [Image Modalities](../run_together_images.md) · [Overview](../run_together.md)
- Step-by-step (manual modules): [Xenium end-to-end tutorial](../../vignettes/pipelines/xenium.md)
