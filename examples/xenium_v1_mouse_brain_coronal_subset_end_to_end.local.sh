#!/bin/bash
set -euo pipefail

##########################################################################################################
# NOTE: This is an end-to-end example of running cartloader on Xenium v1 mouse brain coronal subset data
# To run this script correctly, it is important to specify the absolute paths to the tools
# PLEASE UPDATE the following variables to the absolute paths to the tools
CARTLOADER=/Users/hmkang/code/working/cartloader
PUNKST=${CARTLOADER}/submodules/punkst
CMAP=${CARTLOADER}/assets/fixed_color_map_256.tsv
SPATULA=${CARTLOADER}/submodules/spatula/bin/spatula
TIPPECANOE=${CARTLOADER}/submodules/tippecanoe/tippecanoe
AWS=aws
OUTDIR=./out
##########################################################################################################

URL=https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
ID=xenium-v1-ff-mouse-brain-coronal-subset-ctx-hp

## Download the Xenium output from 10x Genomics
wget ${URL}
## Unzip the downloaded file to the OUTDIR
mkdir -p ${OUTDIR}
unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip -d ${OUTDIR}/outs

## input file directory
INDIR=${OUTDIR}/outs
TSV=${OUTDIR}/tsv/transcripts.unsorted.tsv.gz
FIC=${OUTDIR}/fic
CART=${OUTDIR}/cartl

THREADS=4
JOBS=2
COLNAME=count

## convert the 10x-specific transcript file into generic TSV file compatible with cartloader 
if [ -e "${INDIR}/transcripts.csv.gz" ]; then
	cartloader sge_convert --makefn sge_convert.mk --platform 10x_xenium --in-csv ${INDIR}/transcripts.csv.gz --out-dir ${OUTDIR}/tsv --exclude-feature-regex '^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated)' --sge-visual --spatula ${SPATULA} --n-jobs ${JOBS} --pigz-threads ${THREADS} --csv-colnames-others cell_id z_location
elif [ -e "${INDIR}/transcripts.parquet" ]; then
	cartloader sge_convert --makefn sge_convert.mk --platform 10x_xenium --in-parquet ${INDIR}/transcripts.parquet --out-dir ${OUTDIR}/tsv --exclude-feature-regex '^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated)' --sge-visual --spatula ${SPATULA} --n-jobs ${JOBS} --pigz-threads ${THREADS} --csv-colnames-others cell_id z_location
else
	echo "ERROR: ${INDIR}/transcripts.csv.gz or ${INDIR}/transcripts.parquet files found"
	exit
fi

## write the list of input files to a TSV file
## for multiple input files, each line should be in the format of <sample_name>\t<input_file_path>
echo -e "rep1\t${OUTDIR}/tsv/transcripts.unsorted.tsv.gz" > ${OUTDIR}/tsv/in_list.tsv

## run FICTURE2 with 12um hexagons and 12, 24, 48 factors
LIST=${OUTDIR}/tsv/in_list.tsv
cartloader run_ficture2_multi --in-list ${LIST} --out-dir ${FIC} --width 12 --n-factor 12,24,48 --threads ${THREADS} --n-jobs ${JOBS} --ficture2 ${PUNKST} --cmap-file ${CMAP} --spatula ${SPATULA} --exclude-feature-regex "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated).*" --min-ct-per-unit-hexagon 50 --gzip "pigz -p${THREADS}"

## Write inputn file to import segmented cells, clusters, and cell boundaries from Xenium output
echo -e "rep1\t${INDIR}/analysis/clustering/gene_expression_graphclust/clusters.csv" > ${OUTDIR}/tsv/in_clust.tsv
echo -e "rep1\t${INDIR}/cells.csv.gz" > ${OUTDIR}/tsv/in_xy.tsv
echo -e "rep1\t${INDIR}/cell_boundaries.csv.gz" > ${OUTDIR}/tsv/in_boundaries.tsv

## Use the largest FICTURE model to import cells, clusters, and cell boundaries
MODEL=${FIC}/t12_f48.model.tsv

## Perform LDA-based cell clustering and pixel-level decoding based on segmented cells and boundaries provided by Xenium output
cartloader run_ficture2_multi_cells --all --out-prefix cartloader --out-dir ${FIC} --threads ${THREADS} --n-jobs ${JOBS} --ficture2 ${PUNKST} --cmap-file ${CMAP} --spatula ${SPATULA} --exclude-feature-regex "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated).*" --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --pretrained-model ${MODEL} --gzip "pigz -p${THREADS}"

## Import Xenium Ranger cell clustering and pixel-level decoding
cartloader run_ficture2_multi_cells --all --out-dir ${FIC} --out-prefix xeniumranger --list-cluster ${OUTDIR}/tsv/in_clust.tsv --list-xy ${OUTDIR}/tsv/in_xy.tsv --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --threads 20 --n-jobs 3 --ficture2 ${PUNKST} --cmap-file ${CMAP} --spatula ${SPATULA} --exclude-feature-regex "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated).*" --xy-colname-x x_centroid --xy-colname-y y_centroid --pretrained-model ${MODEL} --gzip "pigz -p${THREADS}"

## Import all 
## We use 'rep1' as generic sample name. Change it to your own sample name that you prefer
SAMPLE=rep1
cartloader run_cartload2 --fic-dir ${FIC}/samples/${SAMPLE} --in-cell-params ${FIC}/samples/${SAMPLE}/ficture.cartloader.params.json ${FIC}/samples/${SAMPLE}/ficture.xeniumranger.params.json --out-dir ${CART}/samples/${SAMPLE} --id ${ID} --n-jobs 5 --threads ${THREADS} --colname-count count --spatula ${SPATULA} --tippecanoe ${TIPPECANOE}

## Import morphology images
FOCUSDIR=${INDIR}/morphology_focus
if [ -d "${FOCUSDIR}" ]; then
	if [ -e "${FOCUSDIR}/morphology_focus_0000.ome.tif" ]; then
		MORPH_NO=0000
		MORPH_TYPE=dapi
		MORPH_COLOR=0F73E6
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		cartloader import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${CART}/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		echo -e "    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles" >> ${CART}/samples/${SAMPLE}/catalog.yaml
	fi
	if [ -e "${FOCUSDIR}/morphology_focus_0001.ome.tif" ]; then
		MORPH_NO=0001
		MORPH_TYPE=boundary
		MORPH_COLOR=F300A5
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		cartloader import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${CART}/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		echo -e "    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles" >> ${CART}/samples/${SAMPLE}/catalog.yaml
	fi
	if [ -e "${FOCUSDIR}/morphology_focus_0002.ome.tif" ]; then
		MORPH_NO=0002
		MORPH_TYPE=rna
		MORPH_COLOR=A4A400
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		cartloader import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${CART}/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		echo -e "    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles" >> ${CART}/samples/${SAMPLE}/catalog.yaml
	fi
	if [ -e "${FOCUSDIR}/morphology_focus_0003.ome.tif" ]; then
		MORPH_NO=0003
		MORPH_TYPE=protein
		MORPH_COLOR=008A00
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		cartloader import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${CART}/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		echo -e "    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles" >> ${CART}/samples/${SAMPLE}/catalog.yaml
	fi
elif [ -e "${INDIR}/morphology_focus.ome.tif" ]; then
	MORPH_TYPE=dapi
	MORPH_COLOR=0F73E6
	TIF=${INDIR}/morphology_focus.ome.tif
	cartloader import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${CART}/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
	echo -e "    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles" >> ${CART}/samples/${SAMPLE}/catalog.yaml
fi

## Upload to S3, Zenodo, or any other storage service
# (grep -E "\." ${CART}/samples/${SAMPLE}/catalog.yaml | perl -lane 'print $F[$#F]' | sort | uniq; echo catalog.yaml;) | xargs -I {} ${AWS} s3 cp ${CART}/samples/${SAMPLE}/{} s3://{BUCKET}/{PREFIX}/{ID}/{}
