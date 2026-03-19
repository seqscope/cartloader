# Quick Start: Run with Docker

This tutorial walks through the `CartLoader` workflow, all packaged inside a prebuilt Docker image with dependencies and input data included.

!!! info "Use Cases"
    This is the fastest and simplest way to try `CartLoader`: no setup, installation, or data download required.

!!! warning "Requirements"
    Users will need to:

    - Set up `Docker` on their system (see [Set Up Docker](#set-up-docker) guide).

---

## Install Docker

If you are new to [Docker](https://www.docker.com/), please refer to the [Docker documentation](https://docs.docker.com/get-started/) for installation and basic usage.

Verify whether Docker is properly set up on your system:

```bash
# Check if Docker is installed and show its version
docker --version

# Test if Docker can successfully run a container
docker run hello-world
```

If these commands fail, [install Docker](https://docs.docker.com/get-docker/) on your system.

---

## Set Up Environment

{%
include-markdown "../../../includes/includemd_vigenettes_setupenv_docker.md"
%}

---

## Prepare Input

### Data Access

!!! warning "Input data already included"
    The [example dataset](https://zenodo.org/records/17953582/files/seqscope_starter.std.tar.gz) is preloaded in the Docker image, so there is no need to download it separately.

    If needed, it is also available on Zenodo: [DOI: 10.5281/zenodo.15701393](https://doi.org/10.5281/zenodo.15701393)

### File Format

The input is a mouse hippocampus SGE in a `FICTURE`-compatible format, prepared by [`sge_convert`](../../reference/sge_convert.md) in `CartLoader`.

{%
  include-markdown "../../../includes/includemd_vigenettes_sgeformat.md"
%}

### Define ID and Parameters

```bash
# Number of jobs
n_jobs=10                                # If not specified, the number of jobs defaults to 1.

# Unique identifier for your dataset
DATA_ID="seqscope_hippo"               # change this to reflect your dataset name
PLATFORM="seqscope"                    # platform information

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor values are provided)
```

---

## SGE Format Conversion

The example dataset is already provided in `FICTURE`-compatible SGE format, so this conversion step is not required in this quickstart.

---

## `FICTURE` Analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
  start="<!--parameter-start-->" end="<!--parameter-end-->"
%}

---

## `CartLoader` Asset Packaging

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
  start="<!--parameter-start-->" end="<!--parameter-end-->"
%}

---

## Upload to Data Repository

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md" start="<!--section1-start-->" end="<!--section1-end-->" preserve-includer-indent=false
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload_aws_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md" start="<!--section2-start-->" end="<!--section2-end-->" preserve-includer-indent=false
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload_zenodo_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md" start="<!--section3-start-->" end="<!--section3-end-->" preserve-includer-indent=false
%}

---

## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_seqscope.md"
%}
