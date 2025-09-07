# 🚀 Quick Start: Run with Docker

This tutorial walks through the `cartloader` workflow — all packaged inside a prebuilt Docker image with dependencies and input data included.

!!! info "Use Cases"

    This is the fastest and simplest way to try `cartloader` — no setup, installation, or data download required.

!!! warning "Requirements"
    Users will need to:

    * Set up `Docker` on their system (see [Set Up Docker](#setting-up-docker) guide).


## Set Up Docker

If you are new to [Docker](https://www.docker.com/), please refer to the [Docker documentation](https://docs.docker.com/get-started/) for installation and basic usage.

Verify whether Docker is properly set up on your system:
```bash
# Check if Docker is installed and show its version
docker --version

# Test if Docker can successfully run a container
docker run hello-world
```

If these commands fail, [install Docker](https://docs.docker.com/get-docker/) in your system.


## Input Data

The input is an mouse hippocampus SGE in a `FICTURE`-compatible format compatible, prepared by [`sge_convert`](../../reference/sge_convert.md) in `cartloader` .

**File Format**


{%
  include-markdown "../../../includes/includemd_vigenettes_sgeformat.md"
%}

**Data Access**

!!! warning "Input data already included"

    The [example dataset](https://zenodo.org/records/15786632/files/seqscope_starter.std.tar.gz) is preloaded in the Docker image — no need to download separately.

    If needs, it is also available on Zenodo: [DOI: 10.5281/zenodo.15701393](https://doi.org/10.5281/zenodo.15701393)

---------------

## Set Up the Environment

!!! info "Fixed paths in the Docker Image"

    Tools and dependencies have fixed paths in the Docker image (e.g., `/usr/local/bin/pmtiles`), which are used directly in the commands below. Skip specifying them manually.

```bash
# ====
# Replace each placeholder with the actual path on your system.  
# ====

work_dir=/path/to/work/directory        # path to work directory that contains the downloaded input data
cd $work_dir

# Number of jobs
n_jobs=10                               # If not specify, the number of jobs defaults to 1.

# Unique identifier for your dataset
DATA_ID="seqscope_hippo"                # change this to reflect your dataset name
PLATFORM="seqscope"                     # platform information

# LDA parameters
train_width=18                           # define LDA training hexagon width (comma-separated if multiple widths are applied)
n_factor=6,12                            # define number of factors in LDA training (comma-separated if multiple n-factor are applied)
```

## `FICTURE` Analysis

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_run_ficture2.md"
  start="<!--parameter-start-->" end="<!--parameter-end-->"
%}

## `cartloader` Compilation

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_run_cartload2.md"
  start="<!--parameter-start-->" end="<!--parameter-end-->"
%}

## Upload to Data Repository

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md"
  start="<!--section1-start-->" end="<!--section1-end-->"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload_aws_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md"
  start="<!--section2-start-->" end="<!--section2-end-->"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload_zenodo_docker.md"
%}

{%
  include-markdown "../../../includes/includemd_vigenettes_upload.md"
  start="<!--section3-start-->" end="<!--section3-end-->"
%}

---------------

## Output Data

{%
  include-markdown "../../../includes/includemd_vigenettes_output_seqscope.md"
%}

