# Which interface should I use to run the pipeline?

CartLoader offers three levels of control. They produce the same kind of CartoScope assets — they differ in how much the tool orchestrates for you.

| Interface | What it is | Best for |
|-----------|------------|----------|
| **`run_together`** | One command runs the whole pipeline (ingest → FICTURE → cells → packaging → images → optional publish) across platforms and samples, as a single resumable Makefile | Multi-platform or multi-sample runs; a joint model across sections; one-command convenience; batch processing |
| **`run_xenium` / `run_visiumhd`** | Per-platform orchestrators for a single platform and sample, with explicit action flags and per-tool path options | A single Xenium/Visium HD sample when you want fine-grained control over each action and tool path |
| **Individual modules** | Call [`sge_convert`](../reference/sge_convert.md), [`run_ficture2`](../reference/run_ficture2.md), [`run_cartload2`](../reference/run_cartload2.md), `import_*`, `upload_*` yourself | Custom pipelines, unsupported/custom input formats, debugging a single step, tight HPC integration |

## Quick decision rule

- **Several samples, or want a joint model, or just want it to run end-to-end?** → `run_together`. Start with the [single-sample](../vignettes/pipelines/run_together.md) or [multi-sample](../vignettes/pipelines/run_together_multi.md) tutorial.
- **One Xenium/Visium HD sample and you want explicit control per step and per tool?** → [`run_xenium`](../reference/run_xenium.md) / [`run_visiumhd`](../reference/run_visiumhd.md).
- **A platform without a built-in profile, a custom format, or you're debugging one step?** → the individual modules. See [Supported Platforms & Inputs](../reference/supported_platforms.md) and the platform starter tutorials.

!!! info "They interoperate"
    All three call the same underlying modules and produce the same `catalog.yaml` + PMTiles. You can prototype with `run_together`, then drop to modules for a step that needs customization — or vice versa.

See also: [Run with Docker or locally?](./choose_run_mode.md) · [Which platform tutorial?](./choose_platform.md)
