# Where are outputs written, and which files should I keep?

## With `run_together`

Everything goes under `--out-dir`:

```
<out-dir>/
├── run_together.mk            # the generated Makefile
├── run_together.resolved.json # the fully resolved configuration (provenance)
├── mk/                        # per-stage done flags
├── tsv/<id>/                  # converted SGE
├── fic/                       # FICTURE models, decoding, and cell analyses
└── cartl/<dir>/               # catalog.yaml + PMTiles  ← the deployable result
```

The packaged directory `<dir>` is `<out-dir name>-<id>`, shortened to `<id>` when the two are equal. A joint multi-sample run writes one such directory per sample plus `cartl/multi-catalog.yaml`.

**Keep** `cartl/` (that is what CartoScope reads; upload it as a unit) and `run_together.resolved.json` (it records exactly what was run). `tsv/` and `fic/` are needed to rerun or extend the analysis, e.g. to add a model with a different `--n-factor` without repeating SGE conversion; delete them only when you are done.

## With the individual commands

Each command writes into its own `--out-dir`, plus a manifest that records its key artifacts:

- [`sge_convert`](../reference/sge_convert.md): unified SGE files + `sge_assets.json`
- [`run_ficture2`](../reference/run_ficture2.md): models and decoding outputs + `ficture.params.json` ([`run_ficture2_multi`](../reference/run_ficture2_multi.md) writes `ficture.multi.params.json`)
- [`run_cartload2`](../reference/run_cartload2.md): PMTiles layers + `catalog.yaml` (and `ficture_assets.json`)

**Keep** the final manifests and `catalog.yaml`, the PMTiles they reference, and any SGE you plan to reuse. **Remove later**, if needed, intermediate files that the pipeline regenerates, and diagnostic plots or logs you no longer need. Before deleting, confirm your next step does not depend on those files.
