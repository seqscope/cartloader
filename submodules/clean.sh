## build htslib
cd htslib
make clean
cd ..

## build qgenlib
cd qgenlib
rm -rf build/
rm -rf bin/
cd ../

## build spatula
cd spatula
rm -rf build/
rm -rf bin/
cd ../

## build pmpoint
cd pmpoint
rm -rf build/
rm -rf bin/
cd ../

## build tippecanoe
cd tippecanoe
make clean
cd ..

## build punkst
cd punkst
rm -rf build/
rm -rf bin/
cd ../

## build ImageMagick
# cd ImageMagick
# make clean
# cd ..

## download pmtiles CLI
rm -rf ./pmtiles
