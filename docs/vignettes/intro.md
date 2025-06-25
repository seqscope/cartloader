# Getting Started with `cartloader`

We offer a series of tutorials to help users master the `cartloader` toolkit across a range of scenarios — from simple getting-started examples to processing full-scale ST datasets. 

For a more detailed explanation of each step, refer to our [References](../reference/sge_convert.md) Page.

## Quick Start

If you're new to `cartloader`, we recommend beginning with the [Quick Start](./quickstart.md) tutorial. 

This vignette provides a beginner-friendly walkthrough using a small mouse hippocampus dataset to walks through `cartloader` functions.

<div class="grid cards single-left" markdown>

-   <figure markdown="span">
    ![Quick Start](../images/starter_vignettes/seqscope.t18_f12_p18_a6.png){ width="100%" }
    <figcaption class="figure-caption-vigintro-large">Quick Start Tutorial</figcaption>
    </figure>

    [Read](./quickstart.md){ .md-button .md-button--primary .button-tight-small }

</div>

## Getting Started per Platform

These tutorials show how to process data from different ST platforms using small, representative subregions.

For each platform, we extract a subset of SGE representing the mouse **hippocampus** as input. If hippocampal data is not available for a given platform, we select one subregion from the mouse brain to illustrate the workflow.

### Sequencing-based Platforms

<div class="grid cards generic" markdown>

-   <figure markdown="span">
    ![example1](../images/starter_vignettes/seqscope.t18_f12_p18_a6.png)
    <figcaption class="figure-caption-vigintro-large">**SeqScope** Starter Guide</figcaption>
    </figure>

    [Read](./quickstart.md){ .md-button .md-button--primary .button-tight-small }

-   <figure markdown="span">
    ![example2](../images/starter_vignettes/visiumhd.t18_f12_p18_a6.png)
    <figcaption class="figure-caption-vigintro-large">**10X Visium HD** Starter Guide</figcaption>
    </figure>

    [Read](./subregion_tutorials/visiumhd.md){ .md-button .md-button--primary .button-tight-small }

-   <figure markdown="span">
    ![example3](../images/starter_vignettes/stereoseq.t18_f12_p18_a6.png){ width="100%", height="300" }
    <figcaption class="figure-caption-vigintro-large">**StereoSeq** Starter Guide</figcaption>
    </figure>

    [Read](./subregion_tutorials/stereoseq.md){ .md-button .md-button--primary .button-tight-small }

-   <figure markdown="span">
    ![example3](../images/starter_vignettes/pixelseq.t18_f12_p18_a6.png){ width="100%" }
    <figcaption class="figure-caption-vigintro-large">**Pixel-seq** Starter Guide</figcaption>
    </figure>

    [Read](./subregion_tutorials/pixelseq.md){ .md-button .md-button--primary .button-tight-small }

</div>

### Imaging-based Platforms
<div class="grid cards generic" markdown>
-   <figure markdown="span">
    ![example3](../images/starter_vignettes/xenium.t12_f12_p12_a6.png){ width="100%" }
    <figcaption class="figure-caption-vigintro-large">**10x Xenium** Starter Guide</figcaption>
    </figure>

    [Read](./subregion_tutorials/xenium.md){ .md-button .md-button--primary .button-tight-small }

-   <figure markdown="span">
    ![example3](../images/starter_vignettes/merscope.t12_f12_p12_a6.png){ width="100%" }
    <figcaption class="figure-caption-vigintro-large">**Vizgen MERSCOPE** Starter Guide</figcaption>
    </figure>

    [Read](./subregion_tutorials/merscope.md){ .md-button .md-button--primary .button-tight-small }

-   <figure markdown="span">
    ![example3](../images/starter_vignettes/cosmxsmi.t12_f12_p12_a6.png){ width="100%" }
    <figcaption class="figure-caption-vigintro-large">**CosMX SMI** Starter Guide</figcaption>
    </figure>

    [Read](./subregion_tutorials/cosmx_smi.md){ .md-button .md-button--primary .button-tight-small }
</div>

## Real-World Use Cases
These examples demonstrate how to process walk through real-world **whole-sample** ST datasets. Each example uses a whole mouse brain sample for consistency.


## Multi-Sample Batch Analysis

!!! warning 
    This function is still under development

This upcoming section will showcase how to process **multiple samples** in a single batch. It will include tutorials on cross-sample analysis.
