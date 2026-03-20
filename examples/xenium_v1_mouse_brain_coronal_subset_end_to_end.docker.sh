#!/bin/bash
set -euo pipefail

##########################################################################################################
# NOTE: This is an end-to-end example of running cartloader on Xenium v1 mouse brain coronal subset data
# To run this script correctly, it is important to specify the absolute paths to the tools
# PLEASE UPDATE the following variables to the absolute paths to the tools
IMAGE="hyunminkang/cartloader:20260319c"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
REAL_OUTDIR="${SCRIPTDIR}/out"
REAL_INDIR="${SCRIPTDIR}/data"

#######################################################################################################
## Set binary paths
CARTLOADER=/app/cartloader
##########################################################################################################

docker pull ${IMAGE}

URL=https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
ID=xenium-v1-ff-mouse-brain-coronal-subset-ctx-hp

## Download the Xenium output from 10x Genomics
wget ${URL}
## Unzip the downloaded file to the OUTDIR
mkdir -p ${REAL_OUTDIR}
unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip -d ${REAL_INDIR}

## input file directory
INDIR=/data
OUTDIR=/out

THREADS=4
JOBS=2
COLNAME=count

CMD="docker run --rm -v ${REAL_OUTDIR}:${OUTDIR} -v ${REAL_INDIR}:${INDIR} ${IMAGE}"
BASH="docker run --rm -v ${REAL_OUTDIR}:${OUTDIR} -v ${REAL_INDIR}:${INDIR} --entrypoint /bin/bash ${IMAGE}"

## convert the 10x-specific transcript file into generic TSV file compatible with cartloader 
if [ -e "${REAL_INDIR}/transcripts.csv.gz" ]; then
	${CMD} sge_convert --makefn sge_convert.mk --platform 10x_xenium --in-csv ${INDIR}/transcripts.csv.gz --out-dir ${OUTDIR}/tsv --exclude-feature-regex '^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated)' --sge-visual --n-jobs ${JOBS} --pigz-threads ${THREADS} --csv-colnames-others cell_id z_location
elif [ -e "${REAL_INDIR}/transcripts.parquet" ]; then
	${CMD} sge_convert --makefn sge_convert.mk --platform 10x_xenium --in-parquet ${INDIR}/transcripts.parquet --out-dir ${OUTDIR}/tsv --exclude-feature-regex '^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated)' --sge-visual --n-jobs ${JOBS} --pigz-threads ${THREADS} --csv-colnames-others cell_id z_location
else
	echo "ERROR: Neither ${INDIR}/transcripts.csv.gz nor ${INDIR}/transcripts.parquet files found"
	exit
fi

## write the list of input files to a TSV file
## for multiple input files, each line should be in the format of <sample_name>\t<input_file_path>
${BASH} -c "echo -e \"rep1\t${OUTDIR}/tsv/transcripts.unsorted.tsv.gz\" > ${OUTDIR}/tsv/in_list.tsv"

## run FICTURE2 with 12um hexagons and 12, 24, 48 factors
LIST=${REAL_OUTDIR}/tsv/in_list.tsv
${CMD} run_ficture2_multi --in-list ${OUTDIR}/tsv/in_list.tsv --out-dir ${OUTDIR}/fic --width 12 --n-factor 12,24,48 --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated).*" --min-ct-per-unit-hexagon 50 --gzip pigz

## Write inputn file to import segmented cells, clusters, and cell boundaries from Xenium output
${BASH} -c "echo -e \"rep1\t${INDIR}/analysis/clustering/gene_expression_graphclust/clusters.csv\" > ${OUTDIR}/tsv/in_clust.tsv"
${BASH} -c "echo -e \"rep1\t${INDIR}/cells.csv.gz\" > ${OUTDIR}/tsv/in_xy.tsv"
${BASH} -c "echo -e \"rep1\t${INDIR}/cell_boundaries.csv.gz\" > ${OUTDIR}/tsv/in_boundaries.tsv"

## Use the largest FICTURE model to import cells, clusters, and cell boundaries
MODEL=${OUTDIR}/fic/t12_f48.model.tsv

## Perform LDA-based cell clustering and pixel-level decoding based on segmented cells and boundaries provided by Xenium output
${CMD} run_ficture2_multi_cells --all --out-prefix cartloader --out-dir ${OUTDIR}/fic --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated).*" --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --pretrained-model ${MODEL} --gzip pigz

## Import Xenium Ranger cell clustering and pixel-level decoding
${CMD} run_ficture2_multi_cells --all --out-dir ${OUTDIR}/fic --out-prefix xeniumranger --list-cluster ${OUTDIR}/tsv/in_clust.tsv --list-xy ${OUTDIR}/tsv/in_xy.tsv --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --threads 20 --n-jobs 3 --exclude-feature-regex "^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated).*" --xy-colname-x x_centroid --xy-colname-y y_centroid --pretrained-model ${MODEL} --gzip pigz

## Import all 
## We use 'rep1' as generic sample name. Change it to your own sample name that you prefer
SAMPLE=rep1
${CMD} run_cartload2 --fic-dir ${OUTDIR}/fic/samples/${SAMPLE} --in-cell-params ${OUTDIR}/fic/samples/${SAMPLE}/ficture.cartloader.params.json ${OUTDIR}/fic/samples/${SAMPLE}/ficture.xeniumranger.params.json --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --id ${ID} --n-jobs 5 --threads ${THREADS} --colname-count count --use-pmpoint --gzip pigz

## Import morphology images
REAL_FOCUSDIR=${REAL_INDIR}/morphology_focus
FOCUSDIR=${INDIR}/morphology_focus
if [ -d "${REAL_FOCUSDIR}" ]; then
	if [ -e "${REAL_FOCUSDIR}/morphology_focus_0000.ome.tif" ]; then
		MORPH_NO=0000
		MORPH_TYPE=dapi
		MORPH_COLOR=0F73E6
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		${CMD} import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		${BASH} -c "echo -e \"    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles\" >> ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml"
	fi
	if [ -e "${REAL_FOCUSDIR}/morphology_focus_0001.ome.tif" ]; then
		MORPH_NO=0001
		MORPH_TYPE=boundary
		MORPH_COLOR=F300A5
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		${CMD} import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		${BASH} -c "echo -e \"    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles\" >> ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml"
	fi
	if [ -e "${REAL_FOCUSDIR}/morphology_focus_0002.ome.tif" ]; then
		MORPH_NO=0002
		MORPH_TYPE=rna
		MORPH_COLOR=A4A400
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		${CMD} import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		${BASH} -c "echo -e \"    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles\" >> ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml"
	fi
	if [ -e "${REAL_FOCUSDIR}/morphology_focus_0003.ome.tif" ]; then
		MORPH_NO=0003
		MORPH_TYPE=protein
		MORPH_COLOR=008A00
		TIF=${FOCUSDIR}/morphology_focus_${MORPH_NO}.ome.tif
		${CMD} import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
		${BASH} -c "echo -e \"    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles\" >> ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml"
	fi
elif [ -e "${REAL_INDIR}/morphology_focus.ome.tif" ]; then
	MORPH_TYPE=dapi
	MORPH_COLOR=0F73E6
	TIF=${INDIR}/morphology_focus.ome.tif
	${CMD} import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5
	${BASH} -c "echo -e \"    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles\" >> ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml"
fi

## Upload to S3, Zenodo, or any other storage service
# (grep -E "\." ${REAL_OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml | perl -lane 'print $F[$#F]' | sort | uniq; echo catalog.yaml;) | xargs -I {} ${AWS} s3 cp ${REAL_OUTDIR}/cartl/samples/${SAMPLE}/{} s3://{BUCKET}/{PREFIX}/{ID}/{}
