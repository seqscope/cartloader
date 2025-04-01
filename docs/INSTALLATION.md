# Installation Guide

This document provides step-by-step instructions to set up the required environment, install necessary dependencies, and install `cartloader`.

---

## 1. Dependencies

To ensure full functionality of cartloader, the following dependencies must be installed:

### **1.1 System Dependencies**  
Ensure the following command-line tools are available on your system:  
- `gzip`  
- `sort`  
- `bgzip`  
- `tabix`  

### **1.2 External Tools & Utilities**  
The following external tools and utilities are required for handling spatial data and file processing. Some of these are included as submodules within the repository.

**Python & Related Packages:**
- [`python`](https://www.python.org/) (`cartloader` has been tested for compatibility with Python v3.10, and v3.13.1)
- [`parquet-tools`](https://github.com/apache/parquet-mr/tree/master/parquet-tools)
- [`ficture`](https://github.com/seqscope/ficture)

**External Tools** (Included in submodules)
- [`spatula`](https://github.com/seqscope/spatula)
- [`tippecanoe`](https://github.com/mapbox/tippecanoe)
- [`magick`](https://imagemagick.org/)
- [`go-pmtiles`](https://github.com/protomaps/go-pmtiles)

**Geospatial Data Handling:**
- [`gdal`](https://gdal.org/)

**Cloud & CLI Tools:**
- [`aws-cli`](https://aws.amazon.com/cli/)

---

## 2. Setting Up the Environment using Conda

We recommended to use Conda to manage dependencies efficiently and avoid conflicts. 

### 2.1 Installing Conda
If Conda is not installed, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

Here is an example to install Miniconda3 on Linux. 

```bash
env_dir=/path/to/your/directory/hosting/tools/
cd $env_dir

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2.2 Creating a Conda Environment
Create a dedicated Conda environment for `cartloader`:
```bash
conda_env=cartloader_env        # replace cartloader_env with the name of your conda environment
python_version=3.13.1           # replace 3.13.1 with the version you prefer

conda create -n cartloader_env python=$python_version
conda activate cartloader_env
```

### 2.3 Installing Dependencies with Conda
Once inside the environment, install the required dependencies:
```bash
conda install -c conda-forge gdal aws-cli imagemagick parquet-tools
```

---

### 3. Installing Cartloader

Clone the repository and install it:

```bash
cd $env_dir
git clone git@github.com:seqscope/cartloader.git
cd cartloader
pip install -e ./                                   # (TODO) Install in editable mode for developmentartloader

# Install its requirements:
pip install -r ./requirements.txt
```

---

### 4. Installing submodules

### 4.1 Installing spatula and its dependencies
Install spatula with its dependencies from the submodules directory:

```bash
cd ${env_dir}/cartloader/submodules
bash -x build.sh
cd ..
```

### 4.2 Installing FICTURE
Install FICTURE and its requirements:
```bash
cd ${env_dir}/cartloader/submodules/ficture
pip install -e ./

# Install additional requirements:
pip install -r requirements.txt
```

### 4.3 Installing tippecanoe

```bash
cd ${env_dir}/cartloader/submodules/tippecanoe
make -j
make install  # Requires root access. 
              # Alternatively, run `make install PREFIX=${env_dir}/cartloader/submodules/tippecanoe/` to specify a custom installation directory.
```

### 4.4 Installing go-pmtiles
An easy way to install go-pmtiles is to download a release from [the official website](https://github.com/protomaps/go-pmtiles/releases) and decompress it. 
This will return a `pmtiles` bin file ready for use.

### 4.5 (Optional) Installing ImageMagic
If you have already installed ImageMagic when [setting conda environment](#23-installing-dependencies-with-conda), skip this step.

```bash
cd ${env_dir}/cartloader/submodules/ImageMagick
./configure # Alternatively, run `./configure --prefix=${env_dir}/cartloader/submodules/ImageMagick`.
make 
make install 
```

---

## 5. Verifying the Installation

To confirm the package is installed correctly, run:
```bash
python -c "import cartloader; print('cartloader installed successfully!')"
```

---
