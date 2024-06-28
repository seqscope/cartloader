## build htslib
cd htslib
autoreconf -i
./configure
make
cd ..

## build qgenlib
cd qgenlib
mkdir build
cd build
cmake ..
make
cd ../../

## build spatula
cd spatula
mkdir build
cd build
cmake ..
make
cd ../../

## build tippecanoe
cd tippecanoe
make
cd ..
