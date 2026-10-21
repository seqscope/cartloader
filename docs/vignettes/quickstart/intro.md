# Get Started

The recommended way to run `CartLoader` is [`run_together`](../../reference/run_together.md): a single command that takes a sample from raw platform output to deployable [CartoScope](https://cartoscope.github.io/) assets: SGE conversion, `FICTURE` spatial factors, cell analyses, background images, and PMTiles packaging. It plans every step as one resumable Makefile, so rerunning the same command after an interruption picks up where it stopped.

Both quick starts below run two examples:

1. **A small mouse hippocampus dataset (Seq-Scope)**, finishes in minutes. In Docker it needs no download at all.
2. **A Xenium human kidney section (GSE264334)**, which exercises the whole pipeline, including cell segmentation imports and the DAPI image.

## Choose how to run

<div class="grid cards forref two-up quickstart" markdown>

-   <ul>
      <li>
        <h3> 🚀 Run with Docker</h3>
        <p>Run <code>cartloader</code> inside the prebuilt image <code>hyunminkang/cartloader</code>.</p>
        <ul>
          <li>✅ No installation required</li>
          <li>📦 All tools and dependencies included</li>
          <li>🧪 Starter dataset preloaded</li>
        </ul>
        <p>
    </ul>
    [Read](./run_in_docker.md){ .md-button .md-button--primary .button-tight-small }

-   <ul>
      <li>
        <h3> 🔧 Run Locally</h3>
        <p>Install <code>cartloader</code> and its tools on your system.</p>
        <ul>
          <li>💻 Use the latest code from GitHub</li>
          <li>⚙️ Full flexibility</li>
          <li>📂 Easy to customize and develop</li>
        </ul>
        <p>
      </li>
    </ul>
    [Read](./run_locally.md){ .md-button .md-button--primary .button-tight-small }

</div>

!!! question "Not sure which to choose?"
    See [Docker or local?](../../faq/choose_run_mode.md). In short: use Docker to try `CartLoader` or to run a released version reproducibly; install locally to follow the latest code or modify it.

## After the quick start

- **Your own data**: the [platform tutorials](../intro.md#platform-tutorials) show the `run_together` command for each supported platform, and [Specifying Inputs](../../reference/run_together_inputs.md) covers every way to describe a sample.
- **Several samples together**: [Multiple Samples with run_together](../pipelines/run_together_multi.md) trains one shared model across samples.
- **Publishing and annotation**: `run_together --s3-upload` and `--anno` / `--anno-deep` (see [Publishing](../../reference/run_together.md#publishing)).
- **Individual commands**: the [Step-by-step](../step_by_step.md) section runs each stage (`sge_convert`, `run_ficture2`, `run_cartload2`, ...) separately for full control.
