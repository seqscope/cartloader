Compute spatial factors using `punkst` (FICTURE2 mode). See more details in [Reference page](../docs/reference/run_ficture2.md).

<!-- ```bash
cartloader run_ficture2 \
    --makefn run_ficture2.mk \                          # (optional) file name of the output make file
    --main \                                            # run all five steps in `run_ficture2`
    --in-transcript ./sge/transcripts.unsorted.tsv.gz \ # path to input transcript-level SGE file
    --in-feature ./sge/feature.clean.tsv.gz \           # (optional) path to input feature file
    --in-minmax ./sge/coordinate_minmax.tsv \           # (optional) path to input minmax file
    --cmap-file ${cmap} \                               # (optional) path to input color map file
    --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \        # regex pattern to exclude features (removing mitochondrial and predicted genes in the example analysis)
    --out-dir ./ficture2 \                              # path to output directory
    --width ${train_width} \                            # LDA training hexagon width (comma-separated if multiple widths are applied)
    --n-factor ${n_factor} \                            # number of factors in LDA training (comma-separated if multiple n-factor are applied)
    --spatula ${spatula} \                              # (optional) path to the spatula binary
    --ficture2 ${punkst} \                              # (optional) path to the punkst directory
    --n-jobs 10  \                                      # (optional) number of parallel jobs 
    --threads 10                                        # (optional) number of threads per job
``` -->


```bash
cartloader run_ficture2 \
  --makefn run_ficture2.mk \
  --main \
  --in-transcript ./sge/transcripts.unsorted.tsv.gz \
  --in-feature ./sge/feature.clean.tsv.gz \
  --in-minmax ./sge/coordinate_minmax.tsv \
  --cmap-file ${cmap} \
  --exclude-feature-regex '^(mt-.*$|Gm\d+$)' \
  --out-dir ./ficture2 \
  --width ${train_width} \
  --n-factor ${n_factor} \
  --spatula ${spatula} \
  --ficture2 ${punkst} \
  --n-jobs 10 \
  --threads 10
```

| Parameter                   | Required  | Type                          | Description                                                                                     |
|-----------------------------|-----------|-------------------------------|-------------------------------------------------------------------------------------------------|
| `--makefn`                  |           | string                        | File name for the generated Makefile                                                            |
| `--main`                    | required  | flag                          | Run all five steps in `run_ficture2`                                                            |
| `--in-transcript`           | required  | string                        | Path to input transcript-level SGE file                                                         |
| `--in-feature`              |           | string                        | Path to input feature file                                                                      |
| `--in-minmax`               |           | string                        | Path to coordinate min/max file                                                                 |
| `--cmap-file`               |           | string                        | Path to color map file                                                                          |
| `--exclude-feature-regex`   |           | regex                         | Pattern to exclude features (e.g., mitochondrial or predicted genes)                            |
| `--out-dir`                 | required  | string                        | Path to output directory                                                                        |
| `--width`                   | required  | int or comma-separated list   | LDA training hexagon width(s)                                                                   |
| `--n-factor`                | required  | int or comma-separated list   | Number of LDA factors                                                                           |
| `--spatula`                 |           | string                        | Path to the `spatula` binary                                                                      |
| `--ficture2`                |           | string                        | Path to the `punkst` directory                                                                    |
| `--n-jobs`                  |           | int                           | Number of parallel jobs                                                                         |
| `--threads`                 |           | int                           | Number of threads per job                                                                       |