#!/bin/bash
set -euo pipefail

##########################################################################################################
# NOTE: This is an end-to-end example of running cartloader (installed locally) on a public Xenium human
# kidney dataset from patients with IgA nephropathy (Bull et al. 2024; GEO GSE264334), reformatted into
# the standard Xenium output layout (transcripts.csv.gz, cells.csv.gz, cell_boundaries.csv.gz,
# analysis/clustering/..., morphology.ome.tif).
#
# `cartloader` must be on your PATH (see the installation guide); its bundled tools are found under
# its submodules/ directory.
#
# The input is downloaded once (about 1.5 GB) into ${INDIR}; delete or move it to download again.
# To use a copy you already have, set INDIR to a directory holding those files before running.
OUTDIR="${OUTDIR:-./out}"
INDIR="${INDIR:-${OUTDIR}/outs}"
##########################################################################################################

URL=https://cartostore.s3.amazonaws.com/data/batch%3D2026_05/xenium-geo-public-dataset-collection/xenium-human-kidney-igan-bull2024-20260516/gse264334-reformatted.tar.gz
ID=xenium-human-kidney-igan-bull2024

## Download and unpack the input (skipped when the transcripts are already in INDIR)
mkdir -p ${OUTDIR} ${INDIR}
if [ ! -e "${INDIR}/transcripts.csv.gz" ]; then
    wget -O ${INDIR}/gse264334-reformatted.tar.gz ${URL}
    tar xzf ${INDIR}/gse264334-reformatted.tar.gz -C ${INDIR} --strip-components=1
    rm ${INDIR}/gse264334-reformatted.tar.gz
fi

THREADS=4
JOBS=2

CMD="cartloader"

REGEX_STR="^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated|System|Gm[0-9]|MT-|mt-|Rps|Rpl|NCS-|NCP-)"
XTRA_COLS="cell_id z_location overlaps_nucleus"
BIN_COUNT=50

WIDTH=12
N_FACTOR="12,24,48"
LARGEST_N_FACTOR=$(echo ${N_FACTOR} | tr ',' '\n' | sort -nr | head -n 1)

## convert the 10x-specific transcript file into a generic TSV file compatible with cartloader
${CMD} sge_convert --makefn sge_convert.mk --platform 10x_xenium --in-csv ${INDIR}/transcripts.csv.gz --out-dir ${OUTDIR}/tsv --exclude-feature-regex "${REGEX_STR}" --sge-visual --n-jobs ${JOBS} --pigz-threads ${THREADS} --csv-colnames-others ${XTRA_COLS} --gzip pigz

## write the list of input files to a TSV file
## for multiple input files, each line should be in the format of <sample_name>\t<input_file_path>
echo -e "rep1\t${OUTDIR}/tsv/transcripts.unsorted.tsv.gz" > ${OUTDIR}/tsv/in_list.tsv

## run FICTURE2 with 12um hexagons and 12, 24, 48 factors
LIST=${OUTDIR}/tsv/in_list.tsv
${CMD} run_ficture2_multi --in-list ${LIST} --out-dir ${OUTDIR}/fic --width ${WIDTH} --n-factor ${N_FACTOR} --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "${REGEX_STR}" --min-ct-per-unit-hexagon 50 --gzip pigz --single-molecule

## Write input files to import segmented cells, clusters, and cell boundaries from the Xenium output
echo -e "rep1\t${INDIR}/analysis/clustering/gene_expression_graphclust/clusters.csv" > ${OUTDIR}/tsv/in_clust.tsv
echo -e "rep1\t${INDIR}/cells.csv.gz" > ${OUTDIR}/tsv/in_xy.tsv
echo -e "rep1\t${INDIR}/cell_boundaries.csv.gz" > ${OUTDIR}/tsv/in_boundaries.tsv

## Use the largest FICTURE model to import cells, clusters, and cell boundaries
MODEL=${OUTDIR}/fic/t${WIDTH}_f${LARGEST_N_FACTOR}.model.tsv

## Perform LDA-based cell clustering and pixel-level decoding based on segmented cells and boundaries provided by the Xenium output
${CMD} run_ficture2_multi_cells --all --out-prefix cartloader --out-dir ${OUTDIR}/fic --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "${REGEX_STR}" --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --pretrained-model ${MODEL} --gzip pigz

## Import Xenium Ranger cell clustering and pixel-level decoding
${CMD} run_ficture2_multi_cells --all --out-dir ${OUTDIR}/fic --out-prefix xeniumranger --list-cluster ${OUTDIR}/tsv/in_clust.tsv --list-xy ${OUTDIR}/tsv/in_xy.tsv --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "${REGEX_STR}" --xy-colname-x x_centroid --xy-colname-y y_centroid --pretrained-model ${MODEL} --gzip pigz

## We use 'rep1' as generic sample name. Change it to your own sample name that you prefer
SAMPLE=rep1
${CMD} run_cartload2 --fic-dir ${OUTDIR}/fic/samples/${SAMPLE} --in-cell-params ${OUTDIR}/fic/samples/${SAMPLE}/ficture.cartloader.params.json ${OUTDIR}/fic/samples/${SAMPLE}/ficture.xeniumranger.params.json --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --id ${ID} --n-jobs ${JOBS} --threads ${THREADS} --colname-count count --pmtiles-pmpoint --gzip pigz --bin-count ${BIN_COUNT}

## Import the DAPI morphology image (a single multi-page OME-TIFF in this dataset)
if [ -e "${INDIR}/morphology.ome.tif" ]; then
    MORPH_TYPE=dapi
    MORPH_COLOR=0F73E6
    TIF=${INDIR}/morphology.ome.tif
    ${CMD} import_image --ome2png --png2pmtiles --georeference --in-img ${TIF} --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --img-id ${MORPH_TYPE} --upper-thres-quantile 0.95 --level 0 --colorize ${MORPH_COLOR} --transparent-below 5 --use-middle-page --high-memory
    echo -e "    ${MORPH_TYPE}: ${MORPH_TYPE}.pmtiles" >> ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml
fi

## Upload to S3, Zenodo, or any other storage service
# (grep -E "\." ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml | perl -lane 'print $F[$#F]' | sort | uniq; echo catalog.yaml;) | xargs -I {} aws s3 cp ${OUTDIR}/cartl/samples/${SAMPLE}/{} s3://{BUCKET}/{PREFIX}/{ID}/{}
