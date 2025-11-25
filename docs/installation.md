# Installation Guide

This document walks through the environment setup and installation steps required for `CartLoader`.

---

## 1. Dependencies

Before installing `CartLoader`, make sure the tools below are present on your system.

### **1.1 Required System Utilities**

Confirm that these command-line programs are installed:

- `gzip`
- `sort`
- `bgzip`
- `tabix`

### **1.2 External Tools and Utilities**

These packages support spatial data handling and file conversion. Several are bundled as git submodules.

**Python & Related Tooling**

- [`python`](https://www.python.org/) (verified with versions *3.10* and *3.13.1*)
- [`parquet-tools`](https://github.com/apache/parquet-mr/tree/master/parquet-tools)

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

We recommend isolating the project in a `conda` environment to avoid dependency conflicts.

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
conda create -n ENV_NAME python=3.13.1   # choose the environment name and Python version that suits your workflow
conda activate ENV_NAME
```

### 2.3 Install Core Dependencies

Install the core dependencies once the environment is active:

```bash
conda install -c conda-forge gdal aws-cli imagemagick parquet-tools r-base
```

---

## 3. Installing `CartLoader`

Clone the repository and install required Python and R packages:

```bash
cd $env_dir
git clone git@github.com:seqscope/cartloader.git
cd cartloader

# Install python requirements:
pip install -r ./requirements.txt

# Install R packages:
Rscript ./install_r_packages.R

pip install -e ./
```

---

## 4. Initializing Submodules
```bash
cd $env_dir/cartloader
git submodule update --init --recursive
```

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

!!! info "What are `FICTURE` and `punkst`?"

    [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2) is a segmentation-free method that infers latent spatial factors—coherent spatial patterns of gene activity—that correspond to underlying transcriptional programs or tissue structures. These factors can then be projected back to the pixel level. Although `FICTURE` is built on a Latent Dirichlet Allocation (LDA) framework by default, it is also compatible with clustering outputs from external tools like `Seurat` for pixel-level projection.

    The [`punkst`](https://github.com/Yichen-Si/punkst) toolkit is a streamlined implementation of [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2), designed for improved computational efficiency and scalability while producing results equivalent to the original.

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

Run the following command to verify that the package loads correctly:

```bash
python -c "import cartloader; print('cartloader installed successfully!')"
```
