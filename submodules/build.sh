#!/bin/bash
set -Eeuo pipefail

## build htslib
cd htslib
autoreconf -i
./configure
make -j 4
cd ..

## build qgenlib
cd qgenlib
mkdir -p build
cd build
cmake ..
make -j 4
cd ../../

## build spatula
cd spatula
mkdir -p build
cd build
cmake ..
make -j 4
cd ../../

## build pmpoint
cd pmpoint
mkdir -p build
cd build
cmake ..
make -j 4
cd ../../

## build tippecanoe
cd tippecanoe
make -j 4
cd ..

## build punkst
cd punkst
mkdir -p build
cd build
cmake ..
make -j 4
cd ../../

## build ImageMagick --- SKIP : build only optionally ---
# cd ImageMagick
# ./configure --without-lqr
# make -j 4
# cd ..

## ficture1 build is also skipped because it will be deprecated

## download pmtiles CLI
bash install-pmtiles.sh --dest ./pmtiles
