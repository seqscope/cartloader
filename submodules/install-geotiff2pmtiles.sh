#!/usr/bin/env bash
#
# install-geotiff2pmtiles.sh
#
# Downloads the geotiff2pmtiles CLI binary from the pspoerri/geotiff2pmtiles
# GitHub releases. Supports Linux and macOS (Darwin) on x86_64 / arm64.
# The upstream project does not publish Windows binaries.
#
# Usage:
#   ./install-geotiff2pmtiles.sh                    # install latest to ./geotiff2pmtiles/
#   ./install-geotiff2pmtiles.sh ./mydir            # install latest to ./mydir/
#   ./install-geotiff2pmtiles.sh --version v0.19    # install specific version
#   ./install-geotiff2pmtiles.sh --dest /usr/local/bin
#
# Environment overrides:
#   GEOTIFF2PMTILES_VERSION   - version tag (e.g. "v0.19")
#   GEOTIFF2PMTILES_DEST      - destination directory (default: <script dir>/geotiff2pmtiles)
#   GEOTIFF2PMTILES_OS        - force OS  (linux | darwin)
#   GEOTIFF2PMTILES_ARCH      - force arch (amd64 | arm64)

set -euo pipefail

REPO="pspoerri/geotiff2pmtiles"
GITHUB_API_BASE="https://api.github.com/repos/${REPO}/releases"
BINARY_NAME="geotiff2pmtiles"

# ── script directory (used for default destination) ────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── defaults ────────────────────────────────────────────────────────────────
VERSION="${GEOTIFF2PMTILES_VERSION:-}"
DEST="${GEOTIFF2PMTILES_DEST:-}"

# ── parse CLI args ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --version|-v) VERSION="$2"; shift 2 ;;
    --dest|-d)    DEST="$2";    shift 2 ;;
    --help|-h)
      sed -n '3,21p' "$0"
      exit 0 ;;
    -*)
      echo "Unknown option: $1" >&2; exit 1 ;;
    *)
      # Treat positional argument as destination directory
      if [[ -z "$DEST" ]]; then
        DEST="$1"
      else
        echo "Unknown argument: $1" >&2; exit 1
      fi
      shift ;;
  esac
done

# Default destination if not specified: <script dir>/geotiff2pmtiles
DEST="${DEST:-${SCRIPT_DIR}/geotiff2pmtiles}"

# ── helper: HTTP GET (works with curl or wget) ─────────────────────────────
http_get() {
  local url="$1"
  if command -v curl &>/dev/null; then
    curl -fsSL "$url"
  elif command -v wget &>/dev/null; then
    wget -qO- "$url"
  else
    echo "Error: curl or wget is required." >&2
    exit 1
  fi
}

http_download() {
  local url="$1" dest="$2"
  if command -v curl &>/dev/null; then
    curl -fSL --progress-bar -o "$dest" "$url"
  elif command -v wget &>/dev/null; then
    wget -q --show-progress -O "$dest" "$url"
  fi
}

# ── detect OS (asset names use lowercase: linux | darwin) ──────────────────
detect_os() {
  if [[ -n "${GEOTIFF2PMTILES_OS:-}" ]]; then echo "$GEOTIFF2PMTILES_OS"; return; fi
  local uname_os
  uname_os="$(uname -s)"
  case "$uname_os" in
    Linux*)    echo "linux"  ;;
    Darwin*)   echo "darwin" ;;
    *)
      echo "Error: unsupported OS '$uname_os'. Upstream releases support only linux and darwin." >&2
      echo "Set GEOTIFF2PMTILES_OS to linux or darwin to override." >&2
      exit 1 ;;
  esac
}

# ── detect CPU architecture (asset names use: amd64 | arm64) ───────────────
detect_arch() {
  if [[ -n "${GEOTIFF2PMTILES_ARCH:-}" ]]; then echo "$GEOTIFF2PMTILES_ARCH"; return; fi
  local uname_arch
  uname_arch="$(uname -m)"
  case "$uname_arch" in
    x86_64|amd64)    echo "amd64" ;;
    aarch64|arm64)   echo "arm64" ;;
    *)
      echo "Error: unsupported architecture '$uname_arch'. Set GEOTIFF2PMTILES_ARCH to amd64 or arm64." >&2
      exit 1 ;;
  esac
}

OS="$(detect_os)"
ARCH="$(detect_arch)"
echo "Detected platform: ${OS}/${ARCH}"

# ── resolve release and find the matching asset via GitHub API ──────────────
if [[ -n "$VERSION" ]]; then
  API_URL="${GITHUB_API_BASE}/tags/${VERSION}"
else
  API_URL="${GITHUB_API_BASE}/latest"
fi

echo "Querying GitHub API for release info..."
RELEASE_JSON="$(http_get "$API_URL")" || {
  echo "Error: failed to fetch release info from GitHub." >&2
  exit 1
}

# Extract tag name
TAG_NAME="$(echo "$RELEASE_JSON" | grep '"tag_name"' | head -1 | cut -d'"' -f4)"
if [[ -z "$TAG_NAME" ]]; then
  echo "Error: could not parse tag_name from GitHub API response." >&2
  exit 1
fi
echo "Release tag: ${TAG_NAME}"

# Extract all browser_download_url values from the release JSON
ASSET_URLS="$(echo "$RELEASE_JSON" | grep '"browser_download_url"' | cut -d'"' -f4)"

if [[ -z "$ASSET_URLS" ]]; then
  echo "Error: release ${TAG_NAME} has no downloadable assets." >&2
  echo "This can happen with tag-only releases. Try specifying a version with --version." >&2
  exit 1
fi

# Find the asset matching: basename == "geotiff2pmtiles-<os>-<arch>"
# Only consider URLs whose final path segment starts with "geotiff2pmtiles-"
# so we don't accidentally pick up the coginfo or pmtransform binaries.
EXPECTED_ASSET="${BINARY_NAME}-${OS}-${ARCH}"
DOWNLOAD_URL=""
while IFS= read -r url; do
  if [[ "$(basename "$url")" == "$EXPECTED_ASSET" ]]; then
    DOWNLOAD_URL="$url"
    break
  fi
done <<< "$ASSET_URLS"

if [[ -z "$DOWNLOAD_URL" ]]; then
  echo "Error: no asset named '${EXPECTED_ASSET}' found in release ${TAG_NAME}." >&2
  echo "" >&2
  echo "Available assets:" >&2
  echo "$ASSET_URLS" | sed 's/^/  /' >&2
  exit 1
fi

ASSET_NAME="$(basename "$DOWNLOAD_URL")"
echo "Asset:        ${ASSET_NAME}"
echo "Download URL: ${DOWNLOAD_URL}"

# ── download ────────────────────────────────────────────────────────────────
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

DOWNLOAD_PATH="${TMPDIR}/${ASSET_NAME}"

echo "Downloading..."
http_download "$DOWNLOAD_URL" "$DOWNLOAD_PATH"

if [[ ! -s "$DOWNLOAD_PATH" ]]; then
  echo "Error: downloaded file is empty: ${DOWNLOAD_PATH}" >&2
  exit 1
fi

# ── install (rename to plain 'geotiff2pmtiles') ────────────────────────────
mkdir -p "$DEST"
mv "$DOWNLOAD_PATH" "${DEST}/${BINARY_NAME}"
chmod +x "${DEST}/${BINARY_NAME}"

INSTALLED_PATH="$(cd "$DEST" && pwd)/${BINARY_NAME}"
echo ""
echo "Installed: ${INSTALLED_PATH}"
echo "Version:   $("${INSTALLED_PATH}" --version 2>/dev/null || "${INSTALLED_PATH}" version 2>/dev/null || echo "${TAG_NAME}")"
echo ""

# ── PATH hint ───────────────────────────────────────────────────────────────
DEST_ABS="$(cd "$DEST" && pwd)"
case ":${PATH}:" in
  *":${DEST_ABS}:"*) ;;
  *)
    echo "Hint: add the binary to your PATH:"
    echo "  export PATH=\"${DEST_ABS}:\$PATH\""
    echo "" ;;
esac

echo "Done."
