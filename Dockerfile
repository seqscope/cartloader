# Use Ubuntu 22.04 as the base ()
FROM ubuntu:22.04

# Use an official Python runtime as a parent image
# FROM python:3.10-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# ===============================
# Install system dependencies
# ===============================

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    cmake \
    pkg-config \
    wget \
    curl \
    gzip \
    pigz \
    perl \
    tabix \
    bc \
    python3 \
    python3-pip \
    python3-dev \
    python-is-python3 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    autoconf \
    automake \
    libtool \
    libdeflate-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgdal-dev \
    libtbb-dev \
    libopencv-dev \
    gdal-bin \
    imagemagick \
    libpng-dev \
    r-base \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Upgrade pip and install aws-cli
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install awscli parquet-tools

# ===============================
# Install cartloader
# ===============================

# Set working directory
WORKDIR /app

# Clone the cartloader repository 
# * submodules: 
#   - Only clone the submodules of punkst, spatula and tippecanoe (skip the rest)
#   - spatula uses docker-dev branch
#   - the following submodules are not installed:
#       - ficture
#       - ImageMagick: (installed in the step apt-get install step)
#       - factor_viz 
#       - go_pmtiles: (download from the release directly)

RUN git clone -b dev --recursive https://github.com/seqscope/cartloader.git 

# Set working directory to the cloned repository
WORKDIR /app/cartloader/submodule

RUN bash -x build.sh

WORKDIR /app/cartloader

# Install Python dependencies
RUN python3 -m pip install --no-cache-dir -r installation/requirements.txt

# Install Python dependencies from ficture 
# * Add this step due to missing python packages
RUN cd assets && \
    wget https://raw.githubusercontent.com/seqscope/ficture/refs/heads/main/requirements.txt --output-document ./ficture_requirements.txt && \
    grep -v '^importlib' ficture_requirements.txt > ficture_requirements.fixed.txt && \
    python3 -m pip install -r ./ficture_requirements.fixed.txt

# Install R dependencies
RUN Rscript installation/install_r_packages.R

# Install cartloader itself
RUN python3 -m pip install -e ./

# ===============================
# Add a test dataset
# ===============================
RUN mkdir -p /app/data && \
    wget https://zenodo.org/records/17953582/files/seqscope_starter.std.tar.gz && \
    tar -xzf seqscope_starter.std.tar.gz -C /app/data && \
    rm seqscope_starter.std.tar.gz

# ===============================
# entrypoint 
# ===============================
# Command to run when starting the container
COPY ./entrypoint.sh /
RUN chmod 755 /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
# Set a default command (fallback if no args passed)
# CMD ["bash"]
CMD []
