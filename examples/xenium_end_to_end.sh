#!/bin/bash
set -euo pipefail

##########################################################################################################
# Generic end-to-end example of running cartloader on arbitrary Xenium output data.
#
# USAGE:
#   ./xenium_end_to_end.sh --id <ID> (--url <URL> | --in-dir <DIR>) [OPTIONS]
#
# REQUIRED ARGUMENTS:
#   --id <str>           Unique identifier for this dataset (used as the cartload --id)
#   --url <url>          URL to the Xenium *_outs.zip file from 10x Genomics (downloaded
#                        and unzipped automatically). Mutually exclusive with --in-dir.
#   --in-dir <dir>       Directory containing already-downloaded/unzipped Xenium output
#                        (e.g. transcripts.csv.gz, cells.csv.gz). Mutually exclusive with --url.
#
# OPTIONAL ARGUMENTS (override defaults, any order):
#   --docker             Run in docker mode (default: local).
#   --width <int>        Hexagon width in um for FICTURE2     (default: 12)
#   --n-factor <list>    Comma-separated list of factors      (default: 12,24,48)
#   --threads <int>      Number of threads per job            (default: 4)
#   --jobs <int>         Number of parallel jobs              (default: 2)
#   --bin-count <int>    Bin count for run_cartload2          (default: 500)
#
# EXAMPLES:
#   ./xenium_end_to_end.sh --id my-dataset \
#       --url https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
#
#   ./xenium_end_to_end.sh --id my-dataset --in-dir /path/to/xenium/outs --docker --threads 8 --jobs 4
#
# NOTE (local): set CARTLOADER below (or have `cartloader` on your PATH).
# NOTE (docker): set IMAGE below to the docker image you want to use.
##########################################################################################################

##########################################################################################################
## Configuration that you may need to update for your environment
IMAGE="hyunminkang/cartloader:latest"   # docker image (only used when SYSTEM=docker)
CARTLOADER="cartloader"                      # cartloader executable/path (only used when SYSTEM=local)
##########################################################################################################

##########################################################################################################
## Parse arguments
usage() {
	echo "Usage: $0 --id <ID> (--url <URL> | --in-dir <DIR>) [--docker] [--width N] [--n-factor LIST] [--threads N] [--jobs N] [--bin-count N]"
}

## Required / mode arguments (no defaults)
ID=""
URL=""
IN_DIR=""
SYSTEM=local

## Optional parameters with default values
WIDTH=12
N_FACTOR=12,24,48
THREADS=4
JOBS=2
BIN_COUNT=500

## Parse named arguments (any order)
while [ "$#" -gt 0 ]; do
	case "$1" in
		--id)        ID=$2;        shift 2 ;;
		--url)       URL=$2;       shift 2 ;;
		--in-dir)    IN_DIR=$2;    shift 2 ;;
		--docker)    SYSTEM=docker; shift 1 ;;
		--width)     WIDTH=$2;     shift 2 ;;
		--n-factor)  N_FACTOR=$2;  shift 2 ;;
		--threads)   THREADS=$2;   shift 2 ;;
		--jobs)      JOBS=$2;      shift 2 ;;
		--bin-count) BIN_COUNT=$2; shift 2 ;;
		-h|--help)   usage; exit 0 ;;
		*)
			echo "ERROR: Unknown argument '$1'"
			usage
			exit 1
			;;
	esac
done

## Validate required arguments
if [ -z "${ID}" ]; then
	echo "ERROR: --id is required"
	usage
	exit 1
fi

if [ -z "${URL}" ] && [ -z "${IN_DIR}" ]; then
	echo "ERROR: one of --url or --in-dir is required"
	usage
	exit 1
fi

if [ -n "${URL}" ] && [ -n "${IN_DIR}" ]; then
	echo "ERROR: --url and --in-dir are mutually exclusive"
	usage
	exit 1
fi

LARGEST_N_FACTOR=$(echo ${N_FACTOR} | tr ',' '\n' | sort -nr | head -n 1)
##########################################################################################################

##########################################################################################################
## Common settings
## Base directory for output:
##   local  -> current working directory (pwd)
##   docker -> the directory containing this script
if [ "${SYSTEM}" == "docker" ]; then
	WORKDIR="$(cd "$(dirname "$0")" && pwd)"
else
	WORKDIR="$(pwd)"
fi
REAL_OUTDIR="${WORKDIR}/out"    # host path holding all cartloader output

## host path holding the Xenium input:
##   --in-dir -> use the provided directory as-is (resolved to an absolute path)
##   --url    -> ${WORKDIR}/data, where the archive will be downloaded and unzipped
if [ -n "${IN_DIR}" ]; then
	if [ ! -d "${IN_DIR}" ]; then
		echo "ERROR: --in-dir '${IN_DIR}' is not a directory"
		exit 1
	fi
	REAL_INDIR="$(cd "${IN_DIR}" && pwd)"
else
	REAL_INDIR="${WORKDIR}/data"
fi

COLNAME=count
REGEX_STR="^(Unassigned|Neg|BLANK|Blank|Intergenic|Deprecated|System|Gm[0-9]|MT-|mt-|Rps|Rpl|NCS-|NCP-)"
XTRA_COLS="cell_id z_location overlaps_nucleus"
SAMPLE=rep1   # generic sample name; change to your preferred sample name
##########################################################################################################

##########################################################################################################
## Set up the execution wrappers and container-visible paths.
## INDIR/OUTDIR are the paths passed to the commands (container paths under docker,
## identical to REAL_INDIR/REAL_OUTDIR under local).
if [ "${SYSTEM}" == "docker" ]; then
	docker pull ${IMAGE}
	INDIR=/data
	OUTDIR=/out
	CMD="docker run --rm -v ${REAL_OUTDIR}:${OUTDIR} -v ${REAL_INDIR}:${INDIR} ${IMAGE}"
	BASH="docker run --rm -v ${REAL_OUTDIR}:${OUTDIR} -v ${REAL_INDIR}:${INDIR} --entrypoint /bin/bash ${IMAGE}"
else
	INDIR=${REAL_INDIR}
	OUTDIR=${REAL_OUTDIR}
	CMD="${CARTLOADER}"
	BASH="bash"
fi
##########################################################################################################

##########################################################################################################
## Download and unzip the Xenium output from 10x Genomics (only when --url is given;
## with --in-dir the input directory is used as-is).
mkdir -p ${REAL_OUTDIR}
if [ -n "${URL}" ]; then
	mkdir -p ${REAL_INDIR}
	ZIPFILE="${WORKDIR}/$(basename ${URL})"
	if [ ! -e "${ZIPFILE}" ]; then
		wget -O ${ZIPFILE} ${URL}
	fi
	## Unzip into REAL_INDIR. Xenium *_outs.zip archives extract their files directly
	## (no wrapping subfolder), so transcripts.* etc. land directly under REAL_INDIR.
	unzip -o ${ZIPFILE} -d ${REAL_INDIR}
fi
##########################################################################################################

## convert the 10x-specific transcript file into generic TSV file compatible with cartloader
if [ -e "${REAL_INDIR}/transcripts.csv.gz" ]; then
	${CMD} sge_convert --makefn sge_convert.mk --platform 10x_xenium --in-csv ${INDIR}/transcripts.csv.gz --out-dir ${OUTDIR}/tsv --exclude-feature-regex "${REGEX_STR}" --sge-visual --n-jobs ${JOBS} --pigz-threads ${THREADS} --csv-colnames-others ${XTRA_COLS} --gzip pigz
elif [ -e "${REAL_INDIR}/transcripts.parquet" ]; then
	${CMD} sge_convert --makefn sge_convert.mk --platform 10x_xenium --in-parquet ${INDIR}/transcripts.parquet --out-dir ${OUTDIR}/tsv --exclude-feature-regex "${REGEX_STR}" --sge-visual --n-jobs ${JOBS} --pigz-threads ${THREADS} --csv-colnames-others ${XTRA_COLS} --gzip pigz
else
	echo "ERROR: Neither ${REAL_INDIR}/transcripts.csv.gz nor ${REAL_INDIR}/transcripts.parquet files found"
	exit 1
fi

## write the list of input files to a TSV file
## for multiple input files, each line should be in the format of <sample_name>\t<input_file_path>
${BASH} -c "echo -e \"${SAMPLE}\t${OUTDIR}/tsv/transcripts.unsorted.tsv.gz\" > ${OUTDIR}/tsv/in_list.tsv"

## run FICTURE2 with the requested hexagon width and factors
LIST=${OUTDIR}/tsv/in_list.tsv
${CMD} run_ficture2_multi --in-list ${LIST} --out-dir ${OUTDIR}/fic --width ${WIDTH} --n-factor ${N_FACTOR} --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "${REGEX_STR}" --min-ct-per-unit-hexagon 50 --gzip pigz --single-molecule

## Write input files to import segmented cells, clusters, and cell boundaries from Xenium output
${BASH} -c "echo -e \"${SAMPLE}\t${INDIR}/analysis/clustering/gene_expression_graphclust/clusters.csv\" > ${OUTDIR}/tsv/in_clust.tsv"
${BASH} -c "echo -e \"${SAMPLE}\t${INDIR}/cells.csv.gz\" > ${OUTDIR}/tsv/in_xy.tsv"
${BASH} -c "echo -e \"${SAMPLE}\t${INDIR}/cell_boundaries.csv.gz\" > ${OUTDIR}/tsv/in_boundaries.tsv"

## Use the largest FICTURE model to import cells, clusters, and cell boundaries
MODEL=${OUTDIR}/fic/t${WIDTH}_f${LARGEST_N_FACTOR}.model.tsv

## Perform LDA-based cell clustering and pixel-level decoding based on segmented cells and boundaries provided by Xenium output
${CMD} run_ficture2_multi_cells --all --out-prefix cartloader --out-dir ${OUTDIR}/fic --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "${REGEX_STR}" --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --pretrained-model ${MODEL} --gzip pigz

## Import Xenium Ranger cell clustering and pixel-level decoding
${CMD} run_ficture2_multi_cells --all --out-dir ${OUTDIR}/fic --out-prefix xeniumranger --list-cluster ${OUTDIR}/tsv/in_clust.tsv --list-xy ${OUTDIR}/tsv/in_xy.tsv --list-boundaries ${OUTDIR}/tsv/in_boundaries.tsv --threads ${THREADS} --n-jobs ${JOBS} --exclude-feature-regex "${REGEX_STR}" --xy-colname-x x_centroid --xy-colname-y y_centroid --pretrained-model ${MODEL} --gzip pigz

## Load everything into cartloader tiles
${CMD} run_cartload2 --fic-dir ${OUTDIR}/fic/samples/${SAMPLE} --in-cell-params ${OUTDIR}/fic/samples/${SAMPLE}/ficture.cartloader.params.json ${OUTDIR}/fic/samples/${SAMPLE}/ficture.xeniumranger.params.json --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --id ${ID} --n-jobs ${JOBS} --threads ${THREADS} --colname-count ${COLNAME} --use-pmpoint --gzip pigz --bin-count ${BIN_COUNT}

##########################################################################################################
## Import morphology images
## import_morph <real_tif> <container_tif> <img_id> <hex_color> [extra import_image args...]
import_morph() {
	local real_tif=$1
	local cont_tif=$2
	local img_id=$3
	local color=$4
	shift 4
	if [ -e "${real_tif}" ]; then
		${CMD} import_image --ome2png --png2pmtiles --georeference --in-img ${cont_tif} --out-dir ${OUTDIR}/cartl/samples/${SAMPLE} --img-id ${img_id} --upper-thres-quantile 0.95 --level 0 --colorize ${color} --transparent-below 5 "$@"
		${BASH} -c "echo -e \"    ${img_id}: ${img_id}.pmtiles\" >> ${OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml"
	fi
}

REAL_FOCUSDIR=${REAL_INDIR}/morphology_focus
FOCUSDIR=${INDIR}/morphology_focus
if [ -d "${REAL_FOCUSDIR}" ]; then
	import_morph "${REAL_FOCUSDIR}/morphology_focus_0000.ome.tif" "${FOCUSDIR}/morphology_focus_0000.ome.tif" dapi     0F73E6
	import_morph "${REAL_FOCUSDIR}/morphology_focus_0001.ome.tif" "${FOCUSDIR}/morphology_focus_0001.ome.tif" boundary F300A5
	import_morph "${REAL_FOCUSDIR}/morphology_focus_0002.ome.tif" "${FOCUSDIR}/morphology_focus_0002.ome.tif" rna      A4A400
	import_morph "${REAL_FOCUSDIR}/morphology_focus_0003.ome.tif" "${FOCUSDIR}/morphology_focus_0003.ome.tif" protein  008A00
elif [ -e "${REAL_INDIR}/morphology_focus.ome.tif" ]; then
	import_morph "${REAL_INDIR}/morphology_focus.ome.tif" "${INDIR}/morphology_focus.ome.tif" dapi 0F73E6
elif [ -e "${REAL_INDIR}/morphology.ome.tif" ]; then
	import_morph "${REAL_INDIR}/morphology.ome.tif" "${INDIR}/morphology.ome.tif" dapi 0F73E6 --use-middle-page --high-memory
fi
##########################################################################################################

## Upload to S3, Zenodo, or any other storage service
# (grep -E "\." ${REAL_OUTDIR}/cartl/samples/${SAMPLE}/catalog.yaml | perl -lane 'print $F[$#F]' | sort | uniq; echo catalog.yaml;) | xargs -I {} aws s3 cp ${REAL_OUTDIR}/cartl/samples/${SAMPLE}/{} s3://{BUCKET}/{PREFIX}/${ID}/{}
