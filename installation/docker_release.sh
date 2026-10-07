#!/bin/bash
# Build, test and (optionally) push the cartloader Docker image for the commit checked out
# here. Run it on a Linux machine (e.g. an EC2 instance) from a fresh clone, inside tmux
# or screen, since a full run takes hours:
#
#   git clone https://github.com/seqscope/cartloader.git && cd cartloader
#   git checkout <branch|tag|commit>      # optional; the default branch otherwise
#   bash installation/docker_release.sh --version 20261007-1432           # build and test
#   bash installation/docker_release.sh --version 20261007-1432 --push    # ... and push
#
# One-time setup on a fresh Amazon Linux 2023 instance:
#   sudo yum install -y docker git tmux awscli-2
#   sudo usermod -a -G docker ec2-user && sudo service docker start
#   newgrp docker                         # or log out and back in
#   docker login                          # only needed for --push
#
# Stages, stopping at the first failure:
#   1. build    The image clones exactly this checkout's commit (which must be on GitHub).
#   2. smoke    installation/docker_smoke_test.sh in the image: dependencies, bundled
#               binaries, and every `cartloader <command>` imports.
#   3. test     `cartloader run_together` on the test dataset, inside the new image
#               (default: the GSE264334 Xenium kidney data, fetched once from S3).
#   4. check    Every file referenced by the output catalog exists and is non-empty.
#   5. compare  The same test with the baseline image (default: the current :latest) and a
#               report of the differences. For review only; it never fails the release.
#   6. push     <repo>:<version> and <repo>:latest. Only with --push, and only if 1-4 passed.
#
# Under --work-dir, runs/<version>/ holds report.md, compare.md, logs/ (one log per stage)
# and out/ (the test outputs). data/ and baseline/ are kept and reused by later runs.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"

usage() {
	cat <<EOF
Usage: $0 --version <tag> [options]
  --version <tag>         Version tag for the image, e.g. 20261007-1432 (required)
  --push                  Push <repo>:<version> and <repo>:latest if all checks pass
  --repo <name>           Image repository (default: hyunminkang/cartloader)
  --skip-build            Reuse <repo>:<version> built by an earlier run of this script
  --baseline-image <img>  Image to compare against (default: <repo>:latest)
  --no-baseline           Skip the baseline run and comparison
  --work-dir <dir>        Working directory (default: \$HOME/cartloader-release)
  --threads <n>           run_together --threads for the test (default: 4)
  --n-jobs <n>            run_together --n-jobs for the test (default: 4)
  --bin-count <n>         run_together --bin-count for the test (default: 50)
  --test-data <s3-uri>    .tar.gz that extracts into a directory of the same name, used as
                          run_together --in-dir (default: the GSE264334 Xenium kidney data)
  --test-platform <name>  run_together --platform for the test data (default: 10x_xenium)
EOF
}

VERSION=""
REPO="hyunminkang/cartloader"
PUSH=0
SKIP_BUILD=0
BASELINE_IMAGE=""
NO_BASELINE=0
WORK="${HOME}/cartloader-release"
THREADS=4
JOBS=4
BIN_COUNT=50
# Public (AWS Open Data) copy of GSE264334 in Xenium folder layout; extracts into gse264334-reformatted/
TEST_DATA="s3://cartostore/data/batch=2026_05/xenium-geo-public-dataset-collection/xenium-human-kidney-igan-bull2024-20260516/gse264334-reformatted.tar.gz"
TEST_PLATFORM="10x_xenium"

while [ "$#" -gt 0 ]; do
	case "$1" in
		--version)        VERSION=$2;        shift 2 ;;
		--push)           PUSH=1;            shift 1 ;;
		--repo)           REPO=$2;           shift 2 ;;
		--skip-build)     SKIP_BUILD=1;      shift 1 ;;
		--baseline-image) BASELINE_IMAGE=$2; shift 2 ;;
		--no-baseline)    NO_BASELINE=1;     shift 1 ;;
		--work-dir)       WORK=$2;           shift 2 ;;
		--threads)        THREADS=$2;        shift 2 ;;
		--n-jobs)         JOBS=$2;           shift 2 ;;
		--bin-count)      BIN_COUNT=$2;      shift 2 ;;
		--test-data)      TEST_DATA=$2;      shift 2 ;;
		--test-platform)  TEST_PLATFORM=$2;  shift 2 ;;
		-h|--help)        usage; exit 0 ;;
		*) echo "ERROR: Unknown argument '$1'"; usage; exit 1 ;;
	esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "${VERSION}" ] || { usage; die "--version is required"; }
IMG="${REPO}:${VERSION}"
BASELINE_IMAGE=${BASELINE_IMAGE:-${REPO}:latest}

##########################################################################################
## Preconditions: everything that would otherwise fail hours into the run
for tool in docker git tar aws; do
	pkg=${tool}; [ "${tool}" = "aws" ] && pkg=awscli-2
	command -v "${tool}" > /dev/null || die "'${tool}' is not installed (sudo yum install -y ${pkg})"
done
docker info > /dev/null 2>&1 || die "cannot reach the docker daemon (is it running, and are you in the docker group?)"

cd "${REPO_DIR}"
SHA=$(git rev-parse HEAD)
if ! git diff --quiet HEAD --; then
	# Only the build reads the working tree (Dockerfile, entrypoint.sh); the image itself
	# clones ${SHA} from GitHub, so a dirty tree would not match the image.
	[ "${SKIP_BUILD}" -eq 1 ] || die "tracked files have uncommitted changes; commit and push them, or check out a clean commit"
	echo "WARNING: tracked files have uncommitted changes (ignored with --skip-build)"
fi
git fetch -q origin || echo "WARNING: git fetch failed; checking against the remote refs from the clone"
if [ -z "$(git branch -r --contains "${SHA}")" ] && [ -z "$(git tag --contains "${SHA}")" ]; then
	die "commit ${SHA} is not on GitHub, but the image clones it from there; push it first"
fi

if [ "${PUSH}" -eq 1 ]; then
	grep -qs "index.docker.io" "${HOME}/.docker/config.json" || die "--push needs registry credentials; run 'docker login' first"
	if docker manifest inspect "${IMG}" > /dev/null 2>&1; then
		die "${IMG} already exists on the registry; choose another --version"
	fi
fi

avail_gb() { df -Pk "$1" | awk 'NR==2 { print int($4 / 1024 / 1024) }'; }
mkdir -p "${WORK}"
WORK="$(cd "${WORK}" && pwd)"
DOCKER_ROOT=$(docker info -f '{{.DockerRootDir}}' 2>/dev/null || echo /var/lib/docker)
for dir in "${WORK}" "${DOCKER_ROOT}"; do
	gb=$(avail_gb "${dir}" 2>/dev/null || echo "")
	if [ -n "${gb}" ] && [ "${gb}" -lt 10 ]; then
		echo "WARNING: only ${gb} GB free under ${dir}; the build and tests may run out of disk"
	fi
done

##########################################################################################
## Run directory. A previous run of the same version is moved aside rather than reused:
## the test's Makefiles would skip finished steps and test nothing. (mv also works on the
## root-owned files the containers leave behind, which rm would not.)
RUN_DIR="${WORK}/runs/${VERSION}"
if [ -e "${RUN_DIR}" ]; then
	mv "${RUN_DIR}" "${RUN_DIR}.$(date +%Y%m%d-%H%M%S)"
fi
mkdir -p "${RUN_DIR}/logs"
TEST_NAME=$(basename "${TEST_DATA}" .tar.gz)   # the archive extracts into this directory
DATA_DIR="${WORK}/data/${TEST_NAME}"
REPORT="${RUN_DIR}/report.md"
T0=${SECONDS}

cat > "${REPORT}" <<EOF
# cartloader image ${IMG}

- commit: ${SHA} ($(git log -1 --format='%cd, %s' --date=short))
- started: $(date '+%Y-%m-%d %H:%M:%S %Z') on $(uname -n), $(nproc) CPUs
- test: run_together --platform ${TEST_PLATFORM} on ${TEST_NAME} (--threads ${THREADS} --n-jobs ${JOBS} --bin-count ${BIN_COUNT})

| stage | result | time |
|---|---|---|
EOF

finish() {
	local status=$1
	{
		echo
		echo "Total time: $(( (SECONDS - T0) / 60 )) min."
		if [ "${status}" -ne 0 ]; then
			echo "**Release checks FAILED; nothing was pushed.**"
		elif [ "${PUSH}" -eq 1 ]; then
			echo "Pushed ${IMG} and ${REPO}:latest."
		else
			echo "All checks passed; nothing was pushed. To publish this image:"
			echo
			echo '```'
			echo "docker tag ${IMG} ${REPO}:latest"
			echo "docker push ${IMG}"
			echo "docker push ${REPO}:latest"
			echo '```'
		fi
	} >> "${REPORT}"
	echo
	cat "${REPORT}"
	echo
	echo "Report: ${REPORT}"
	if [ -f "${RUN_DIR}/compare.md" ]; then
		echo "Comparison with ${BASELINE_IMAGE}: ${RUN_DIR}/compare.md"
	fi
	exit "${status}"
}

## run_stage <name> <fatal: 1|0> <description> <command...>
## Runs the command with its output in logs/<name>.log and records the result. The command
## runs in a subshell with errexit, outside any `if`, so a failing step inside a stage
## function stops that stage (errexit is ignored inside an if condition).
run_stage() {
	local name=$1 fatal=$2 desc=$3; shift 3
	local log="${RUN_DIR}/logs/${name}.log" start=${SECONDS} rc
	echo "[$(date +%H:%M:%S)] ${desc} (log: ${log})"
	set +e
	( set -e; "$@" ) > "${log}" 2>&1
	rc=$?
	set -e
	local mins=$(( (SECONDS - start) / 60 ))
	if [ "${rc}" -eq 0 ]; then
		echo "    passed (${mins} min)"
		echo "| ${name} | passed | ${mins} min |" >> "${REPORT}"
		return 0
	fi
	echo "    FAILED (${mins} min); last lines of the log:"
	tail -n 40 "${log}" | sed 's/^/    /'
	if [ "${fatal}" -eq 1 ]; then
		echo "| ${name} | **FAILED** | ${mins} min |" >> "${REPORT}"
		finish 1
	fi
	echo "| ${name} | failed (not required) | ${mins} min |" >> "${REPORT}"
	return 1
}

##########################################################################################
## Stages

verify_image_commit() {
	local inner
	inner=$(docker run --rm --entrypoint git "${IMG}" -C /app/cartloader rev-parse HEAD)
	if [ "${inner}" != "${SHA}" ]; then
		echo "${IMG} contains commit ${inner}, expected ${SHA}"
		return 1
	fi
	echo "${IMG} contains commit ${SHA}"
}

build_image() {
	docker build --build-arg CARTLOADER_REF="${SHA}" -t "${IMG}" "${REPO_DIR}"
	verify_image_commit
	docker image inspect -f 'image size: {{.Size}} bytes' "${IMG}"
}

smoke_test() {
	docker run --rm -v "${REPO_DIR}/installation:/release:ro" --entrypoint bash "${IMG}" /release/docker_smoke_test.sh
}

fetch_test_data() {
	# The marker sits next to the data, not in it, where run_together would see it.
	if [ -f "${DATA_DIR}.extracted" ]; then
		echo "reusing ${DATA_DIR}"
		return 0
	fi
	local tgz="${WORK}/data/$(basename "${TEST_DATA}")"
	mkdir -p "${WORK}/data"
	# Public data: --no-sign-request needs no AWS credentials on the instance.
	aws s3 cp --no-sign-request --only-show-errors "${TEST_DATA}" "${tgz}.part"
	mv "${tgz}.part" "${tgz}"
	rm -rf "${DATA_DIR}"
	tar -xzf "${tgz}" -C "${WORK}/data"
	rm -f "${tgz}"
	if [ ! -d "${DATA_DIR}" ]; then
		echo "$(basename "${tgz}") did not extract into ${DATA_DIR}"
		return 1
	fi
	touch "${DATA_DIR}.extracted"
}

## run_test <image> <out dir>: run_together on the test data inside the image. The data is
## mounted read-only, so the new and baseline runs share it without affecting each other.
run_test() {
	mkdir -p "$2"
	docker run --rm -v "${DATA_DIR}:/data:ro" -v "$2:/out" "$1" run_together \
		--platform "${TEST_PLATFORM}" --in-dir /data --out-dir /out --id "${TEST_NAME}" \
		--threads "${THREADS}" --n-jobs "${JOBS}" --bin-count "${BIN_COUNT}"
}

check_outputs() {
	docker run --rm -v "${REPO_DIR}/installation:/release:ro" -v "${RUN_DIR}/out:/out:ro" \
		--entrypoint python3 "${IMG}" /release/compare_test_outputs.py check /out
}

BASE_DIR=""   # set in the main flow once the baseline image is pulled
run_baseline() {
	echo "outputs in ${BASE_DIR}"
	if [ -f "${BASE_DIR}/.complete" ]; then
		echo "reusing the finished baseline run"
		return 0
	fi
	if [ -e "${BASE_DIR}" ]; then
		mv "${BASE_DIR}" "${BASE_DIR}.incomplete.$(date +%Y%m%d-%H%M%S)"
	fi
	run_test "${BASELINE_IMAGE}" "${BASE_DIR}/out"
	touch "${BASE_DIR}/.complete"
}

compare_outputs() {
	echo "/new = ${RUN_DIR}/out, /baseline = ${BASE_DIR}/out (${BASELINE_IMAGE})"
	docker run --rm -v "${REPO_DIR}/installation:/release:ro" -v "${RUN_DIR}/out:/new:ro" \
		-v "${BASE_DIR}/out:/baseline:ro" -v "${RUN_DIR}:/report" \
		--entrypoint python3 "${IMG}" /release/compare_test_outputs.py compare /new /baseline --out /report/compare.md
}

push_image() {
	docker tag "${IMG}" "${REPO}:latest"
	docker push "${IMG}"
	docker push "${REPO}:latest"
}

##########################################################################################
## Main
run_stage data 1 "Fetching test data into ${DATA_DIR}" fetch_test_data
if [ "${SKIP_BUILD}" -eq 1 ]; then
	docker image inspect "${IMG}" > /dev/null 2>&1 || die "--skip-build: ${IMG} does not exist locally"
	run_stage build 1 "Checking that the existing ${IMG} contains commit ${SHA}" verify_image_commit
else
	run_stage build 1 "Building ${IMG} from commit ${SHA}" build_image
fi
run_stage smoke 1 "Smoke test inside ${IMG}" smoke_test
run_stage test 1 "run_together on ${TEST_NAME} with ${IMG}" run_test "${IMG}" "${RUN_DIR}/out"
run_stage check 1 "Checking the files referenced by the output catalog" check_outputs

if [ "${NO_BASELINE}" -eq 0 ] && run_stage baseline-pull 0 "Pulling the baseline image ${BASELINE_IMAGE}" docker pull "${BASELINE_IMAGE}"; then
	# Keyed by image ID, so a rerun against the same baseline reuses its finished run.
	BASE_ID=$(docker image inspect -f '{{.Id}}' "${BASELINE_IMAGE}" | cut -d: -f2 | cut -c1-12)
	BASE_DIR="${WORK}/baseline/${TEST_NAME}-${BASE_ID}"
	if run_stage baseline 0 "Same test with the baseline image ${BASELINE_IMAGE} (${BASE_ID})" run_baseline; then
		run_stage compare 0 "Comparing outputs with the baseline" compare_outputs || true
	fi
fi

if [ "${PUSH}" -eq 1 ]; then
	run_stage push 1 "Pushing ${IMG} and ${REPO}:latest" push_image
fi
finish 0
