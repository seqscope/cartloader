#!/bin/bash
# Build, test and (optionally) push the cartloader Docker image for the commit checked out
# here. Run it on a Linux machine (e.g. an EC2 instance) from a fresh clone, inside tmux
# or screen, since a full run takes hours:
#
#   git clone https://github.com/seqscope/cartloader.git && cd cartloader
#   git checkout <branch|tag|commit>      # optional; the default branch otherwise
#   bash installation/docker_release.sh --version 20261007a           # build and test
#   bash installation/docker_release.sh --version 20261007a --push    # ... and push
#
# One-time setup on a fresh Amazon Linux 2023 instance:
#   sudo yum install -y docker git wget unzip tmux
#   sudo usermod -a -G docker ec2-user && sudo service docker start
#   newgrp docker                         # or log out and back in
#   docker login                          # only needed for --push
#
# Stages, stopping at the first failure:
#   1. build    The image clones exactly this checkout's commit (which must be on GitHub).
#   2. smoke    installation/docker_smoke_test.sh in the image: dependencies, bundled
#               binaries, and every `cartloader <command>` imports.
#   3. e2e      examples/xenium_end_to_end.sh --docker against the new image.
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
  --version <tag>         Version tag for the image, e.g. 20261007a (required)
  --push                  Push <repo>:<version> and <repo>:latest if all checks pass
  --repo <name>           Image repository (default: hyunminkang/cartloader)
  --skip-build            Reuse <repo>:<version> built by an earlier run of this script
  --baseline-image <img>  Image to compare against (default: <repo>:latest)
  --no-baseline           Skip the baseline run and comparison
  --work-dir <dir>        Working directory (default: \$HOME/cartloader-release)
  --threads <n>           Threads per job for the test run (default: 4)
  --jobs <n>              Parallel jobs for the test run (default: 2)
  --test-url <url>        Xenium *_outs.zip used for the test (default: mouse brain subset)
  --test-id <id>          Dataset ID for the test (default: matches the default URL)
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
JOBS=2
TEST_URL="https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip"
TEST_ID="xenium-v1-ff-mouse-brain-coronal-subset-ctx-hp"

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
		--jobs)           JOBS=$2;           shift 2 ;;
		--test-url)       TEST_URL=$2;       shift 2 ;;
		--test-id)        TEST_ID=$2;        shift 2 ;;
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
for tool in docker git wget unzip; do
	command -v "${tool}" > /dev/null || die "'${tool}' is not installed (sudo yum install -y ${tool})"
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
DATA_DIR="${WORK}/data/${TEST_ID}"
REPORT="${RUN_DIR}/report.md"
T0=${SECONDS}

cat > "${REPORT}" <<EOF
# cartloader image ${IMG}

- commit: ${SHA} ($(git log -1 --format='%cd, %s' --date=short))
- started: $(date '+%Y-%m-%d %H:%M:%S %Z') on $(uname -n), $(nproc) CPUs
- test: ${TEST_ID} (--threads ${THREADS} --jobs ${JOBS})

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
	if [ -f "${DATA_DIR}/.unzipped" ]; then
		echo "reusing ${DATA_DIR}"
		return 0
	fi
	local zip="${WORK}/data/$(basename "${TEST_URL}")"
	mkdir -p "${WORK}/data"
	if [ ! -f "${zip}" ]; then
		wget -nv -O "${zip}.part" "${TEST_URL}"
		mv "${zip}.part" "${zip}"
	fi
	rm -rf "${DATA_DIR}"
	mkdir -p "${DATA_DIR}"
	unzip -q -o "${zip}" -d "${DATA_DIR}"
	touch "${DATA_DIR}/.unzipped"
}

## run_e2e <image> <work dir>: the documented Xenium example, in docker mode
run_e2e() {
	bash "${REPO_DIR}/examples/xenium_end_to_end.sh" --id "${TEST_ID}" --in-dir "${DATA_DIR}" \
		--docker --image "$1" --no-pull --work-dir "$2" --threads "${THREADS}" --jobs "${JOBS}"
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
	run_e2e "${BASELINE_IMAGE}" "${BASE_DIR}"
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
run_stage e2e 1 "End-to-end test (${TEST_ID}) with ${IMG}" run_e2e "${IMG}" "${RUN_DIR}"
run_stage check 1 "Checking the files referenced by the output catalog" check_outputs

if [ "${NO_BASELINE}" -eq 0 ] && run_stage baseline-pull 0 "Pulling the baseline image ${BASELINE_IMAGE}" docker pull "${BASELINE_IMAGE}"; then
	# Keyed by image ID, so a rerun against the same baseline reuses its finished run.
	BASE_ID=$(docker image inspect -f '{{.Id}}' "${BASELINE_IMAGE}" | cut -d: -f2 | cut -c1-12)
	BASE_DIR="${WORK}/baseline/${TEST_ID}-${BASE_ID}"
	if run_stage baseline 0 "Same test with the baseline image ${BASELINE_IMAGE} (${BASE_ID})" run_baseline; then
		run_stage compare 0 "Comparing outputs with the baseline" compare_outputs || true
	fi
fi

if [ "${PUSH}" -eq 1 ]; then
	run_stage push 1 "Pushing ${IMG} and ${REPO}:latest" push_image
fi
finish 0
