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

RUN git clone https://github.com/seqscope/cartloader.git 

# Set working directory to the cloned repository
WORKDIR /app/cartloader

# sync up submodules 
RUN git submodule update --init submodules/punkst submodules/spatula submodules/tippecanoe

# Install Python dependencies
RUN python3 -m pip install --no-cache-dir -r requirements.txt

# Install Python dependencies from ficture 
# * Add this step due to missing python packages
RUN cd assets && \
    wget https://raw.githubusercontent.com/seqscope/ficture/refs/heads/main/requirements.txt --output-document ./ficture_requirements.txt && \
    grep -v '^importlib' ficture_requirements.txt > ficture_requirements.fixed.txt && \
    python3 -m pip install -r ./ficture_requirements.fixed.txt

# Install cartloader itself
RUN python3 -m pip install -e ./

# ===============================
# Install submodules
# ===============================
RUN cd submodules/spatula && \
    git checkout docker-dev && \
    git submodule update --init --recursive submodules/htslib submodules/qgenlib  

# Build submodule: htslib
RUN cd submodules/spatula/submodules/htslib && \
    autoreconf -i && \
    ./configure && \
    make -j$(nproc)

# Build submodule: qgenlib
RUN cd submodules/spatula/submodules/qgenlib && \
    mkdir -p build && cd build && \
    cmake .. && \
    make -j$(nproc)

# Build submodule: spatula
RUN cd submodules/spatula && \
    mkdir -p build && cd build && \
    cmake .. && \
    make -j$(nproc)

# Build submodule: tippecanoe
# Build 
RUN cd submodules/tippecanoe && \
    make -j && \
    make install

# Build submodule: punkst
# * libtbb-dev and libopencv-dev are installed at the step of installing system dependencies
# * add this sed cmd to update the markerselection.hpp file to include <optional> to avoid the error of `/app/cartloader/submodules/punkst/src/markerselection.hpp:173:14: error: 'optional' is not a member of 'std'`
RUN cd submodules/punkst && \
    git checkout main && \
    git pull origin main && \
    git submodule update --init

RUN sed -i '/#include <tbb\/global_control.h>/a #include <optional>' submodules/punkst/src/markerselection.hpp
RUN cd submodules/punkst && \
    mkdir -p build && cd build && \
    cmake .. && \
    cmake --build .  --parallel 2

# Build tools: go-pmtiles
RUN wget https://github.com/protomaps/go-pmtiles/releases/download/v1.28.0/go-pmtiles_1.28.0_Linux_x86_64.tar.gz && \
    mkdir -p /opt/go-pmtiles && \
    tar -zxvf go-pmtiles_1.28.0_Linux_x86_64.tar.gz --one-top-level=/opt/go-pmtiles && \
    mv /opt/go-pmtiles/pmtiles /usr/local/bin/ && \
    rm -rf go-pmtiles_1.28.0_Linux_x86_64.tar.gz /opt/go-pmtiles


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
