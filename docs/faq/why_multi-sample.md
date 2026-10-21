# Why use Multi‑Sample FICTURE Analysis?

**`CartLoader` supports analyzing ≥2 samples in three ways**:

- **A joint `run_together` run** ([`--samples`](../vignettes/pipelines/run_together_multi.md)): lists the samples in a sheet and trains one shared model across all of them, then packages each sample against it. This is the recommended way; internally it uses `run_ficture2_multi`.
- **Multi‑sample FICTURE analysis with the individual commands** ([`run_ficture2_multi`](../reference/run_ficture2_multi.md)): the same joint model, run step by step.
- **SGE stitch + single‑sample analysis** ([`sge_stitch`](../reference/sge_stitch.md) → [`run_ficture2`](../reference/run_ficture2.md)): stitch multiple SGEs into a single mosaic, then train one model on that mosaic.

For unrelated samples that should **not** share a model, use `run_together --out-root` instead: each sample gets its own independent model.

**What to expect**

- Shared factors/comparability: a joint model learns a cohort‑wide latent basis and returns per‑sample decodes for direct comparison. The stitch approach yields a single model over the merged mosaic; useful when you need a unified coordinate system (e.g., tiling adjacent sections).
- Efficiency and scale: a joint model fits once for the cohort and decodes per sample, avoiding repeated runs and post‑hoc alignment. Stitching can be simpler for mosaics but often increases I/O and memory due to very large merged files.

**Recommendation**

- Prefer a joint `run_together` run (or `run_ficture2_multi`) for most cohorts: clean per‑sample outputs and better computational efficiency. Use stitching when a single shared coordinate frame is required.
- If you choose stitching, plan for higher resource usage (RAM, disk, and I/O). Large mosaics can be slow to generate and train on, and may require substantially more memory and temporary storage than per‑sample runs.
