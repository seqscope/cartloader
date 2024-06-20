#!/bin/bash

# Variables
input=$1
output=$2
cell_ID=$3
Count=$4

# AWK script
awk -v cell_ID="$cell_ID" -v Count="$Count" 'BEGIN { OFS="\t"; print "X", "Y", "gene", cell_ID, Count }
    NR > 1 {
        if ($1 == prevX && $2 == prevY) {
            sumCount += $Count;
        } else {
            if (NR > 2) {
            print prevX, prevY, prevGene, firstCellID, sumCount;
            }
            prevX = $1; prevY = $2; prevGene = $3; firstCellID = $cell_ID; sumCount = $Count;
        }
    }
    END { print prevX, prevY, prevGene, firstCellID, sumCount; }' ${input} | gzip -c > $output
