# Which interface should I use to run the pipeline?

**Use [`run_together`](../reference/run_together.md).** It is the recommended interface: one command runs the whole pipeline (SGE conversion → FICTURE → cell analyses → images → packaging → optional publishing) for any supported platform, for one sample or many, as a single resumable Makefile.

The individual commands remain available when you need control that `run_together` does not expose.

| Interface | What it is | Best for |
|-----------|------------|----------|
| **`run_together`** (recommended) | One command for the whole pipeline, with per-platform defaults, sample sheets, and a JSON config for full control | Almost every run: single samples, cohorts with a joint model, batches of independent samples |
| **Individual commands** | Call [`sge_convert`](../reference/sge_convert.md), [`run_ficture2`](../reference/run_ficture2.md), [`run_cartload2`](../reference/run_cartload2.md), `import_*`, `upload_*` yourself | Custom pipelines, a parameter `run_together` does not expose, debugging one stage, tight HPC integration |
| **`run_xenium` / `run_visiumhd`** (legacy) | Older per-platform orchestrators | Existing scripts only; see [Legacy Orchestrators](../reference/legacy.md) |

## Quick decision rule

- **Running a supported platform, one sample or many?** → `run_together`. Start with [Get Started](../vignettes/quickstart/intro.md), then the [single-sample](../vignettes/pipelines/run_together.md) or [multi-sample](../vignettes/pipelines/run_together_multi.md) guide.
- **A platform without a built-in profile, or a custom format?** → still `run_together`, with `--platform generic` (see the [Generic platform page](../reference/platforms/generic.md)); drop to the individual commands only if that is not enough.
- **Need to customize or debug one stage?** → the individual commands, introduced in [Step-by-step](../vignettes/step_by_step.md).

!!! info "They interoperate"
    `run_together` calls the same commands and produces the same `catalog.yaml` + PMTiles. You can prototype with `run_together`, then rerun a single stage by hand.

See also: [Run with Docker or locally?](./choose_run_mode.md) · [Which platform tutorial?](./choose_platform.md)
