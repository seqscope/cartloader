## build htslib
cd htslib
autoreconf -i
./configure
make -j 4
cd ..

## build qgenlib
cd qgenlib
mkdir build
cd build
cmake ..
make -j 4
cd ../../

## build spatula
cd spatula
mkdir build
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
mkdir build
cd build
cmake ..
make -j 4
cd ../../

## build ImageMagick
cd ImageMagick
./configure --without-lqr
make -j 4
cd ..

## download pmtiles CLI
bash install-pmtiles.sh --dest ./pmtiles
