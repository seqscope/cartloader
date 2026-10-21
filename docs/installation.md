# Installation Guide

This page installs `CartLoader` on your own machine.

!!! tip "No installation needed with Docker"
    The prebuilt image `hyunminkang/cartloader:latest` contains `CartLoader` and every tool it uses. If you only want to run the pipeline, see [Quick Start with Docker](./vignettes/quickstart/run_in_docker.md) instead.

In short:

1. Create an environment for the non-Python tools ([Section 2](#2-setting-up-the-environment-using-conda)).
2. Install `CartLoader` and its Python dependencies from the repository with [`uv`](https://docs.astral.sh/uv/) or `pip` ([Section 3](#3-installing-cartloader)).
3. Build the bundled tools with `build.sh` ([Section 4](#4-initializing-and-building-submodules)).
4. Verify ([Section 5](#5-verifying-the-installation)).

---

## 1. Dependencies

Instructions for installing everything below are in Sections [2](#2-setting-up-the-environment-using-conda), [3](#3-installing-cartloader) and [4](#4-initializing-and-building-submodules).

### **1.1 Required System Utilities**

Confirm that these command-line programs are installed:

- `make`, `gzip`, `pigz`, `sort`, `bc`, `perl`
- `bgzip`, `tabix` (from [htslib](https://www.htslib.org/))

### **1.2 Build Tools**

`build.sh` ([Section 4](#4-initializing-and-building-submodules)) compiles several tools from source, which needs:

- A C/C++ compiler, `cmake`, `autoconf`, `automake`, `libtool`, `pkg-config`, `git`
- Development headers for `zlib`, `bzip2`, `liblzma`, `libcurl`, `openssl`, `libpng`, `libdeflate`, `TBB` and `OpenCV`

On Ubuntu/Debian, the [Dockerfile](https://github.com/seqscope/cartloader/blob/main/Dockerfile) lists the exact `apt` packages. On macOS, install them with Homebrew.

### **1.3 External Tools and Utilities**

**Python & Related Tooling**

- [`python`](https://www.python.org/) 3.10 or newer (verified with versions *3.10* and *3.13*)
- [`uv`](https://docs.astral.sh/uv/) (optional, recommended): a fast drop-in for `pip`

The Python packages `CartLoader` needs, including [`parquet-tools`](https://pypi.org/project/parquet-tools/), are declared in its `pyproject.toml` and installed automatically in [Section 3](#3-installing-cartloader).

**R & Related Packages**

- [`R` from CRAN](https://cran.r-project.org/) (verified with version *4.5.1*). The R packages are installed in [Section 3.2](#32-installing-r-packages).

**Tools Built from the Repository** (git submodules, built by `build.sh`)

- [`punkst`](https://github.com/Yichen-Si/punkst): the latest, more efficient implementation of [FICTURE](https://github.com/seqscope/ficture)
- [`spatula`](https://github.com/seqscope/spatula)
- [`pmpoint`](https://github.com/seqscope/pmpoint)
- [`tippecanoe`](https://github.com/felt/tippecanoe)
- [`go-pmtiles`](https://github.com/protomaps/go-pmtiles) and `geotiff2pmtiles`: prebuilt binaries that `build.sh` downloads

**Other Tools**

- [ImageMagick](https://imagemagick.org/): `convert` (ImageMagick 6) or `magick` (ImageMagick 7); `spatula` uses it for image I/O
- [`gdal`](https://gdal.org/)
- [`aws-cli`](https://aws.amazon.com/cli/), only for uploading to S3

---

## 2. Setting Up the Environment using `conda`

We recommend isolating the project in a `conda` environment to avoid dependency conflicts. `conda` provides the Python interpreter and the non-Python tools (`gdal`, `R`, ...); the Python packages are installed from `pyproject.toml` in [Section 3](#3-installing-cartloader).

!!! tip "Without `conda`"

    If `gdal`, `R` and the other tools in [Section 1](#1-dependencies) already come from your system (e.g. `apt`, Homebrew, or environment modules), you can skip this section and use a plain virtual environment instead; see [Section 3.3](#33-alternative-a-virtual-environment-without-conda).

### 2.1 Installing `conda`

If `conda` is not already available, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

Example installation of `Miniconda3` on Linux:

```bash
env_dir=/path/to/your/directory/hosting/tools/      ## replace with your preferred tools directory
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
conda install -c conda-forge gdal aws-cli imagemagick pigz r-base
```

Install only these tools with `conda`. Leave the Python packages (including `parquet-tools`) to `uv` or `pip` in the next section: having the same Python package installed by both `conda` and `pip` is a common source of conflicts.

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
git clone --recursive https://github.com/seqscope/cartloader.git   # --recursive also fetches the submodules (Section 4)
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

Then build the bundled tools with `build.sh`, the same step the Docker image uses:

```bash
cd $env_dir/cartloader/submodules
bash build.sh
```

It builds `htslib`, `qgenlib`, `spatula`, `pmpoint`, `tippecanoe` and `punkst`, and downloads the `pmtiles` and `geotiff2pmtiles` binaries. ImageMagick is not built; install it with `conda` ([Section 2.3](#23-install-core-dependencies)) or see [Section 4.5](#45-installing-imagemagick).

!!! warning "The tools must stay inside `submodules/`"
    `CartLoader` commands, and `run_together` in particular, look for these tools at fixed locations inside the checkout (for example `submodules/punkst/bin/punkst` and `submodules/pmtiles/pmtiles`), not on your `PATH`. Build them in place as below; a copy installed elsewhere is not found.

If a step of `build.sh` fails, rebuild that tool on its own as follows, then rerun `build.sh` or continue with the remaining tools.

### 4.1 Installing `spatula`

```bash
cd ${env_dir}/cartloader/submodules/spatula
mkdir -p build && cd build
cmake ..
make -j 4
```

`build.sh` builds `htslib` and `qgenlib` first, which `spatula` and `pmpoint` link against; build them first if you start from a clean checkout.

### 4.2 Installing `punkst`

[`punkst`](https://github.com/Yichen-Si/punkst) provides [`FICTURE` (Si et al., Nature Methods 2024)](https://www.nature.com/articles/s41592-024-02415-2). Build it inside the submodule so that `submodules/punkst/bin/punkst` exists:

```bash
cd ${env_dir}/cartloader/submodules/punkst
mkdir -p build && cd build
cmake ..
make -j 4
```

For build prerequisites and troubleshooting, see the [`punkst` installation guide](https://yichen-si.github.io/punkst/install/).

!!! question

    - [What is `FICTURE` and `punkst`?](./faq/ficture.md)

### 4.3 Installing `tippecanoe`

```bash
cd ${env_dir}/cartloader/submodules/tippecanoe
make -j 4
```

`CartLoader` uses the binary built in place (`submodules/tippecanoe/tippecanoe`); `make install` is not needed.

### 4.4 Installing `go-pmtiles` and `geotiff2pmtiles`

Both are downloaded as prebuilt binaries for your OS and architecture:

```bash
cd ${env_dir}/cartloader/submodules
bash install-pmtiles.sh --dest ./pmtiles
bash install-geotiff2pmtiles.sh --dest ./geotiff2pmtiles
```

### 4.5 Installing ImageMagick

Skip this step if ImageMagick was already installed via `conda` in [Section 2.3](#23-install-core-dependencies) or by your system's package manager.

```bash
cd ${env_dir}/cartloader/submodules/ImageMagick
./configure     # Alternatively, run `./configure --prefix=${env_dir}/cartloader/submodules/ImageMagick`.
make
make install
```

---

## 5. Verifying the Installation

Run the following commands to verify `CartLoader` and its dependencies:

```bash
# activate your environment if you did not
conda activate ${ENV_NAME}       # or, for Section 3.3: source .venv/bin/activate

# verify the installation of cartloader: lists the available commands
cartloader

# check the system tools, Python packages and R packages
python ./installation/check_dependencies.py

# check that the bundled tools were built
ls -l submodules/punkst/bin/punkst submodules/spatula/bin/spatula submodules/pmpoint/bin/pmpoint \
      submodules/tippecanoe/tippecanoe submodules/pmtiles/pmtiles submodules/geotiff2pmtiles/geotiff2pmtiles
```

`check_dependencies.py` reads the Python dependencies from `pyproject.toml`. Packages of the optional `[ai]` extra are reported as optional and do not fail the check.

Then try the [local quick start](./vignettes/quickstart/run_locally.md).

---

## 6. Updating `CartLoader`

```bash
cd $env_dir/cartloader
git pull
git submodule update --init --recursive

# pick up any dependency changes in pyproject.toml (include [ai] if you installed it)
uv pip install -e ".[ai]"          # or: python -m pip install -e ".[ai]"

# rebuild the bundled tools for the updated submodule versions
cd submodules && bash build.sh
```

Both steps matter after each update: new code may need a Python package added to `pyproject.toml` (a stale environment fails with an error such as `requires the 'imagecodecs' package`), and the pipeline passes options that only the tool versions pinned by the repository understand (an outdated build fails with an "unknown option" error or a usage message).
