Compute spatial factors using `punkst` (FICTURE2). See more details on the [Reference page](../docs/reference/run_ficture2.md).

```bash
cartloader run_ficture2 \
  --makefn run_ficture2.mk \
  --main \
  --in-transcript ./sge/transcripts.unsorted.tsv.gz \
  --in-feature ./sge/feature.clean.tsv.gz \
  --in-minmax ./sge/coordinate_minmax.tsv \
  --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \
  --out-dir ./ficture2 \
  --width ${train_width} \
  --n-factor ${n_factor} \
  --n-jobs ${n_jobs} \
  --threads ${n_jobs}
```

<!--parameter-start-->
| Parameter                 | Required              | Type                        | Description |
|---------------------------|-----------------------|-----------------------------|-------------|
| `--main`                  | required <sup>1</sup> | flag                        | Run all [five main steps](../docs/reference/run_ficture2.md#actions) |
| `--in-transcript`         | required              | string                      | Path to the transcript-indexed SGE file from `sge_convert` |
| `--out-dir`               | required              | string                      | Output directory |
| `--width`                 | required              | int or comma-separated list | Hexagon width(s) in µm for LDA training |
| `--n-factor`              | required              | int or comma-separated list | Number(s) of LDA factors; one model per width × factor combination |
| `--makefn`                |                       | string                      | File name for the generated Makefile (default: `run_ficture2.mk`) |
| `--in-feature`            |                       | string                      | Path to the per-gene count file from `sge_convert` (default: the feature list made during tiling) |
| `--in-minmax`             |                       | string                      | Path to the coordinate min/max file from `sge_convert` (default: the range computed during tiling) |
| `--exclude-feature-regex` |                       | regex                       | Feature names to leave out of the factor analysis. `punkst` matches the **whole** name, so end each alternative with `.*` or `$` as needed (the example drops mitochondrial `mt-*` and predicted `Gm<digits>` genes) |
| `--cmap-file`             |                       | string                      | Color map TSV (default: `assets/default_color_map.tsv` in the `CartLoader` repository) |
| `--n-jobs`                |                       | int                         | Number of parallel jobs (default: `1`) |
| `--threads`               |                       | int                         | Number of threads per job (default: `4`) |

<sub><sup>1</sup>: `CartLoader` requires **at least one action**: `--tile` (tiling), `--segment` (hexagon segmentation), `--init-lda` (LDA training), `--decode` (pixel-level decoding), `--umap` (UMAP of the factors), or `--main` for all five.</sub>
<!--parameter-end-->
