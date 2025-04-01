#!/bin/bash

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_geotiff> <output_dimensions_file>"
    exit 1
fi

INPUT_GEOTIFF="$1"
OUTPUT_DIMENSIONS_FILE="$2"

# Use gdalinfo to extract dimensions
INFO=$(gdalinfo "$INPUT_GEOTIFF" 2>/dev/null)
if [[ $INFO =~ Size\ is\ ([0-9]+),\ ([0-9]+) ]]; then
    WIDTH=${BASH_REMATCH[1]}
    HEIGHT=${BASH_REMATCH[2]}
    echo -e "WIDTH\t$WIDTH\nHEIGHT\t$HEIGHT" > "$OUTPUT_DIMENSIONS_FILE"
else
    echo "Error: Failed to extract dimensions from $INPUT_GEOTIFF"
    exit 1
fi
