# Why use Multi‑Sample FICTURE Analysis?

**`CartLoader` supports analyzing ≥2 samples in two ways**:

- Multi‑sample FICTURE analysis ([`run_ficture2_multi`](../reference/run_ficture2_multi.md)): Jointly learns spatial factors across all samples and writes per‑sample outputs in one parallelizable run.
- SGE stitch + single‑sample analysis ([`sge_stitch`](../reference/sge_stitch.md) → [`run_ficture2`](../reference/run_ficture2.md)): Stitch multiple SGEs into a single mosaic, then train one model on that mosaic.

**What to expect**

- Shared factors/comparability: `run_ficture2_multi` learns a cohort‑wide latent basis and returns per‑sample decodes for direct comparison. The stitch approach yields a single model over the merged mosaic; useful when you need a unified coordinate system (e.g., tiling adjacent sections).
- Efficiency and scale: `run_ficture2_multi` fits once for the cohort and decodes per sample, avoiding repeated runs and post‑hoc alignment. Stitching can be simpler for mosaics but often increases I/O and memory due to very large merged files.

**Recommendation:**

    - Prefer `run_ficture2_multi` for most cohorts for clean per‑sample outputs and better computational efficiency; use stitching when a single shared coordinate frame is required.
- If you choose stitching, plan for higher resource usage (RAM, disk, and I/O). Large mosaics can be slow to generate and train on, and may require substantially more memory and temporary storage than per‑sample runs.
