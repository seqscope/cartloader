#!/bin/bash

# Check if input file is provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_geotiff>"
    exit 1
fi

# Input GeoTIFF file
GEOTIFF="$1"

# Check if the file exists
if [[ ! -f "$GEOTIFF" ]]; then
    echo "Error: File '$GEOTIFF' not found."
    exit 1
fi

# Extract dimensions using gdalinfo
INFO=$(gdalinfo "$GEOTIFF" 2>&1)
if [[ $INFO =~ Size\ is\ ([0-9]+),\ ([0-9]+) ]]; then
    WIDTH=${BASH_REMATCH[1]}
    HEIGHT=${BASH_REMATCH[2]}
    echo "$WIDTH $HEIGHT"
else
    echo "Error: Failed to extract dimensions from '$GEOTIFF'."
    exit 1
fi

