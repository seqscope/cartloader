# Installation Guide

This document walks through the environment setup and installation steps required for `CartLoader`.

In short: create an environment for the non-Python tools (Section 2), install `CartLoader` and its Python dependencies from the repository with [`uv`](https://docs.astral.sh/uv/) or `pip` (Section 3), build the bundled tools (Section 4), and verify (Section 5).

---

## 1. Dependencies

Below listed all required tools and packages for `CartLoader`. Instruction of how to install those tools and packages are provided in the sections [2](#2-setting-up-the-environment-using-conda), [3](#3-installing-cartloader) and [4](#4-initializing-and-building-submodules).

### **1.1 Required System Utilities**

Confirm that these command-line programs are installed:

- `gzip`
- `sort`
- `bgzip`
- `tabix`
- `bc`
- `perl`

### **1.2 External Tools and Utilities**

These packages support spatial data handling and file conversion. Several are bundled as git submodules.

**Python & Related Tooling**

- [`python`](https://www.python.org/) 3.10 or newer (verified with versions *3.10* and *3.13*)
- [`uv`](https://docs.astral.sh/uv/) (optional, recommended): a fast drop-in for `pip`
- [`parquet-tools`](https://github.com/apache/parquet-mr/tree/master/parquet-tools)

The Python packages that `CartLoader` needs are declared in its `pyproject.toml` and installed automatically in [Section 3](#3-installing-cartloader).

**R & related Packages:**

- [`R` from CRAN](https://cran.r-project.org/)(verified with versions *4.5.1*)
 
**External Tools** (included as submodules)

- [`punkst`](https://github.com/Yichen-Si/punkst) (the latest and more efficient implementation of [FICTURE](https://github.com/seqscope/ficture))
- [`spatula`](https://github.com/seqscope/spatula)
- [`tippecanoe`](https://github.com/mapbox/tippecanoe)
- [`magick`](https://imagemagick.org/)
- [`go-pmtiles`](https://github.com/protomaps/go-pmtiles)

**Geospatial Utilities**

- [`gdal`](https://gdal.org/)

**Cloud & CLI Utilities**

- [`aws-cli`](https://aws.amazon.com/cli/)

---

## 2. Setting Up the Environment using `conda`

We recommend isolating the project in a `conda` environment to avoid dependency conflicts. `conda` provides the Python interpreter and the non-Python tools (`gdal`, `R`, ...); the Python packages are installed from `pyproject.toml` in [Section 3](#3-installing-cartloader).

!!! tip "Without `conda`"

    If `gdal`, `R` and the other tools in [Section 1](#1-dependencies) already come from your system (e.g. `apt`, Homebrew, or environment modules), you can skip this section and use a plain virtual environment instead; see [Section 3.3](#33-alternative-a-virtual-environment-without-conda).

### 2.1 Installing `conda`

If `conda` is not already available, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

Example installation of `Miniconda3` on Linux:

```bash
env_dir=/path/to/your/directory/hosting/tools/      ## replace `/path/to/your/directory/hosting/tools/` with your preferred tools directory
cd $env_dir

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2.2 Creating an Environment

Set up a fresh environment for `CartLoader`:

```bash
ENV_NAME="cart_env"   # define your environment name
PY_VERSION="3.13"     # define your python version (3.10 or newer)

conda create -n ${ENV_NAME} python=${PY_VERSION}
conda activate ${ENV_NAME}
```

### 2.3 Install Core Dependencies

Install the non-Python dependencies once the environment is active:

```bash
conda install -c conda-forge gdal aws-cli imagemagick parquet-tools r-base
```

Install only these tools with `conda`. Leave the Python packages to `uv` or `pip` in the next section: having the same Python package installed by both `conda` and `pip` is a common source of conflicts.

### 2.4 Installing `uv` (optional)

[`uv`](https://docs.astral.sh/uv/) installs Python packages much faster than `pip` and uses the same commands (`uv pip install ...`). Skip this step to use `pip`.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh    # or: pip install uv
```

---

## 3. Installing `CartLoader`

### 3.1 Installing the Python Package

Clone the repository, including its submodules, and install `CartLoader` with its Python dependencies. Run the install inside the activated environment from [Section 2](#2-setting-up-the-environment-using-conda):

```bash
cd $env_dir
git clone --recursive git@github.com:seqscope/cartloader.git   # --recursive also fetches the submodules (Section 4)
cd cartloader
```

=== "uv"

    ```bash
    # uv installs into the active conda environment (it detects $CONDA_PREFIX)
    uv pip install -e .

    # or, with the optional AI annotation packages (see below):
    uv pip install -e ".[ai]"
    ```

=== "pip"

    ```bash
    python -m pip install -e .

    # or, with the optional AI annotation packages (see below):
    python -m pip install -e ".[ai]"
    ```

The optional `[ai]` extra adds the packages for AI annotation of factors and cell clusters: `run_together --anno` / `--anno-deep`, `anno_cartload_folder`, and the `annotate_*` commands. Those commands also need the API key of the provider they call, set as an environment variable (`ANTHROPIC_API_KEY`, `GEMINI_API_KEY`, `OPENAI_API_KEY` or `UMGPT_API_KEY`).

!!! note "Why an editable install (`-e`)"

    `CartLoader` reads its `assets/` and the submodule binaries from the repository checkout, so it must be installed in editable mode, pointing at the checkout, rather than as a regular package. A `git pull` then updates the code without reinstalling; re-run the install command when the dependencies in `pyproject.toml` change (see [Section 6](#6-updating-cartloader)).

### 3.2 Installing R Packages

```bash
Rscript ./installation/install_r_packages.R
```

### 3.3 Alternative: a Virtual Environment without `conda`

If the non-Python tools come from your system instead of `conda`, create a virtual environment in the checkout and install into it:

=== "uv"

    ```bash
    cd $env_dir/cartloader
    uv venv --python 3.13          # creates .venv; uv downloads that Python if it is missing
    source .venv/bin/activate
    uv pip install -e ".[ai]"      # or: uv pip install -e .
    ```

=== "pip"

    ```bash
    cd $env_dir/cartloader
    python3 -m venv .venv          # python3 must be 3.10 or newer
    source .venv/bin/activate
    python -m pip install -e ".[ai]"   # or: python -m pip install -e .
    ```

Then install the R packages as in [Section 3.2](#32-installing-r-packages).

---

## 4. Initializing and Building Submodules

The `--recursive` clone in [Section 3](#3-installing-cartloader) already fetched the submodules. For a checkout cloned without it, fetch them now:

```bash
cd $env_dir/cartloader
git submodule update --init --recursive
```

To build all the bundled tools in one step (the same step the Docker image uses), run `build.sh`. It builds `htslib`, `qgenlib`, `spatula`, `pmpoint`, `tippecanoe` and `punkst`, and downloads the `pmtiles` and `geotiff2pmtiles` binaries; ImageMagick is not built (install it with `conda` in [Section 2.3](#23-install-core-dependencies), or see [Section 4.5](#45-installing-imagemagick)).

```bash
cd $env_dir/cartloader/submodules
bash build.sh
```

If a step fails, or to build a tool on its own, follow the sections below.

### 4.1 Installing `spatula`

Install `spatula` with its dependencies from the submodules directory:

```bash
cd ${env_dir}/cartloader/submodules/spatula

cd submodules
bash -x build.sh
cd ..

## build spatula
mkdir build
cd build
cmake ..
make
```

### 4.2 Installing `punkst`

Install the [`punkst`](https://github.com/Yichen-Si/punkst) toolkit to use [`FICTURE` (Si et al., Nature Methods 2024)](https://www.nature.com/articles/s41592-024-02415-2).

Please follow the [`punkst` installation guide](https://yichen-si.github.io/punkst/install/).

!!! question

    - [What is `FICTURE` and `punkst`?](./faq/ficture.md)

### 4.3 Installing `tippecanoe`

```bash
cd ${env_dir}/cartloader/submodules/tippecanoe
make -j

## Choose one of the following installation options:
# (1) System-wide installation (requires root access):
make install

# (2) Local installation (no root access): specify a custom PREFIX
make install PREFIX=${env_dir}/cartloader/submodules/tippecanoe/  # Replace with your desired installation path
```

### 4.4 Installing `go-pmtiles`
An easy way to install `go-pmtiles` is to download a release from [the official website](https://github.com/protomaps/go-pmtiles/releases) and decompress it. This provides a `pmtiles` binary ready for use.

Here is an example of its installation:

```bash
cd ${env_dir}
wget https://github.com/protomaps/go-pmtiles/releases/download/v1.28.0/go-pmtiles_1.28.0_Linux_x86_64.tar.gz
tar -zxvf go-pmtiles_1.28.0_Linux_x86_64.tar.gz
```

### 4.5 Installing ImageMagick

Skip this step if ImageMagick was already installed via `conda` in [Section 2.3](#23-install-core-dependencies).

```bash
cd ${env_dir}/cartloader/submodules/ImageMagick
./configure     # Alternatively, run `./configure --prefix=${env_dir}/cartloader/submodules/ImageMagick`.
make 
make install 
```

---

## 5. Verifying the Installation

Run the following commands to verify `CartLoader` and all dependencies:

```bash
# activate your environment if you did not
conda activate ${ENV_NAME}       # or, for Section 3.3: source .venv/bin/activate

# verify the installation of cartloader: lists the available commands
cartloader

# verify all dependencies
python ./installation/check_dependencies.py
```

`check_dependencies.py` reads the Python dependencies from `pyproject.toml`. Packages of the optional `[ai]` extra are reported as optional and do not fail the check.

---

## 6. Updating `CartLoader`

```bash
cd $env_dir/cartloader
git pull
git submodule update --init --recursive

# pick up any dependency changes in pyproject.toml (include [ai] if you installed it)
uv pip install -e ".[ai]"          # or: python -m pip install -e ".[ai]"
```

Rebuild the submodules ([Section 4](#4-initializing-and-building-submodules)) when their versions change.
