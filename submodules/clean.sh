## build htslib
cd htslib
make clean
cd ..

## build qgenlib
cd qgenlib
cd build
rm -f CMakeCache.txt
make clean
cd ../../

## build spatula
cd spatula
cd build
rm -f CMakeCache.txt
make clean
cd ../../

## build tippecanoe
cd tippecanoe
make clean
cd ..
