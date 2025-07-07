# Installation Guide

This is an instruction to set up the required environment, and install dependencies and `cartloader`.

---

## 1. Dependencies

To ensure full functionality of `cartloader`, the following dependencies must be installed:

### **1.1 System Dependencies**  

Ensure the following command-line tools are available on your system:  

- `gzip`  
- `sort`  
- `bgzip`  
- `tabix`  

### **1.2 External Tools & Utilities**  
The following external tools and utilities are required for handling spatial data and file processing. Some of these are included as submodules within the repository.

**Python & Related Packages:**

- [`python`](https://www.python.org/) (`cartloader` has been tested with *v3.10* and *v3.13.1*)
- [`parquet-tools`](https://github.com/apache/parquet-mr/tree/master/parquet-tools)

**External Tools** (Included in submodules)

- [`punkst`](https://github.com/Yichen-Si/punkst) (the latest and more efficient implementation of [FICTURE](https://github.com/seqscope/ficture))
- [`spatula`](https://github.com/seqscope/spatula)
- [`tippecanoe`](https://github.com/mapbox/tippecanoe)
- [`magick`](https://imagemagick.org/)
- [`go-pmtiles`](https://github.com/protomaps/go-pmtiles)

**Geospatial Data Handling:**

- [`gdal`](https://gdal.org/)

**Cloud & CLI Tools:**

- [`aws-cli`](https://aws.amazon.com/cli/)

---

## 2. Setting Up the Environment using `conda`

We recommended to use `conda` to manage dependencies efficiently and avoid conflicts.

### 2.1 Installing `conda`

If `conda` is not installed, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

Here is an example to install `Miniconda3` on Linux.

```bash
env_dir=/path/to/your/directory/hosting/tools/      ## replace `/path/to/your/directory/hosting/tools/` by the path to your tool directory
cd $env_dir

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2.2 Creating an Environment

Create a dedicated Conda environment for `cartloader`:

```bash
conda_env=ENV_NAME              # replace ENV_NAME with the name of your conda environment
python_version=3.13.1           # replace 3.13.1 with the version you prefer

conda create -n ENV_NAME python=$python_version
conda activate ENV_NAME
```

### 2.3 Installing Dependencies

Once inside the environment, install the required dependencies:

```bash
conda install -c conda-forge gdal aws-cli imagemagick parquet-tools
```

---

## 3. Installing `cartloader`

Clone the repository and install it:

```bash
cd $env_dir
git clone git@github.com:seqscope/cartloader.git
cd cartloader

pip install -e ./

# Install its requirements:
pip install -r ./requirements.txt
```

---

## 4. Installing Submodules
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

Install [`punkst`](https://github.com/Yichen-Si/punkst) toolkit to use [`FICTURE` (Si et al., Nature Methods 2024)](https://www.nature.com/articles/s41592-024-02415-2).

Please following the [`punkst` installation guide](https://yichen-si.github.io/punkst/install/).

??? "What is `FICTURE` and `punkst`"

    [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2) is a segmentation-free method that infers latent spatial factors—coherent spatial patterns of gene activity—that correspond to underlying transcriptional programs or tissue structures. These factors can then be projected back to the pixel level. Although `FICTURE` is built on a Latent Dirichlet Allocation (LDA) framework by default, it is also compatible with clustering outputs from external tools like `Seurat` for pixel-level projection.

    The [`punkst`](https://github.com/Yichen-Si/punkst) toolkit is a streamlined implementation of the [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2), designed for improved computational efficiency and scalability while producing results equivalent to the original [`FICTURE`](https://www.nature.com/articles/s41592-024-02415-2).

### 4.3 Installing `tippecanoe`

```bash
cd ${env_dir}/cartloader/submodules/tippecanoe
make -j

## Choose one of the following methods to install tippecanoe:
# (1) System-wide installation (requires root access):
make install

# (2) Local installation (no root access): specify a custom PREFIX
make install PREFIX=${env_dir}/cartloader/submodules/tippecanoe/  # Replace with your desired installation path
```

### 4.4 Installing `go-pmtiles`
An easy way to install `go-pmtiles` is to download a release from [the official website](https://github.com/protomaps/go-pmtiles/releases) and decompress it.
This will return a `pmtiles` bin file ready for use.

Here is an example of its installation:

```bash
cd ${env_dir}
wget https://github.com/protomaps/go-pmtiles/releases/download/v1.28.0/go-pmtiles_1.28.0_Linux_x86_64.tar.gz ./
tar -zxvf ./go-pmtiles_1.28.0_Linux_x86_64.tar.gz
```

### 4.5 Installing ImageMagick

If you have already installed ImageMagick when [setting conda environment](#23-installing-dependencies), skip this step.

```bash
cd ${env_dir}/cartloader/submodules/ImageMagick
./configure     # Alternatively, run `./configure --prefix=${env_dir}/cartloader/submodules/ImageMagick`.
make 
make install 
```

---

## 5. Verifying the Installation

To confirm the package is installed correctly, run:

```bash
python -c "import cartloader; print('cartloader installed successfully!')"
```
