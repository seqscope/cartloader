#!/usr/bin/env bash
#
# install-pmtiles.sh
#
# Downloads the pmtiles CLI binary from the go-pmtiles GitHub releases.
# Supports Linux, macOS (Darwin), and Windows (via Git Bash / WSL / MSYS2).
#
# Usage:
#   ./install-pmtiles.sh                    # install latest version to ./bin/
#   ./install-pmtiles.sh ./pmtiles          # install latest to ./pmtiles/
#   ./install-pmtiles.sh --version v1.28.3  # install specific version
#   ./install-pmtiles.sh --dest /usr/local/bin
#
# Environment overrides:
#   PMTILES_VERSION   - version tag (e.g. "v1.28.3")
#   PMTILES_DEST      - destination directory (default: ./bin)
#   PMTILES_OS        - force OS  (Linux | Darwin | Windows)
#   PMTILES_ARCH      - force arch (x86_64 | arm64)

set -euo pipefail

REPO="protomaps/go-pmtiles"
GITHUB_API_BASE="https://api.github.com/repos/${REPO}/releases"

# ── defaults ────────────────────────────────────────────────────────────────
VERSION="${PMTILES_VERSION:-}"
DEST="${PMTILES_DEST:-}"

# ── parse CLI args ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --version|-v) VERSION="$2"; shift 2 ;;
    --dest|-d)    DEST="$2";    shift 2 ;;
    --help|-h)
      sed -n '3,15p' "$0"
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

# Default destination if not specified
DEST="${DEST:-./bin}"

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

# ── detect OS ───────────────────────────────────────────────────────────────
detect_os() {
  if [[ -n "${PMTILES_OS:-}" ]]; then echo "$PMTILES_OS"; return; fi
  local uname_os
  uname_os="$(uname -s)"
  case "$uname_os" in
    Linux*)                    echo "Linux"   ;;
    Darwin*)                   echo "Darwin"  ;;
    CYGWIN*|MINGW*|MSYS*)      echo "Windows" ;;
    *)
      echo "Error: unsupported OS '$uname_os'. Set PMTILES_OS to Linux, Darwin, or Windows." >&2
      exit 1 ;;
  esac
}

# ── detect CPU architecture ─────────────────────────────────────────────────
detect_arch() {
  if [[ -n "${PMTILES_ARCH:-}" ]]; then echo "$PMTILES_ARCH"; return; fi
  local uname_arch
  uname_arch="$(uname -m)"
  case "$uname_arch" in
    x86_64|amd64)    echo "x86_64" ;;
    aarch64|arm64)   echo "arm64"  ;;
    *)
      echo "Error: unsupported architecture '$uname_arch'. Set PMTILES_ARCH to x86_64 or arm64." >&2
      exit 1 ;;
  esac
}

OS="$(detect_os)"
ARCH="$(detect_arch)"
echo "Detected platform: ${OS}/${ARCH}"

# ── resolve release and find the matching asset via GitHub API ──────────────
#
# Instead of guessing the asset filename, we query the API and search
# the actual asset list.  This handles naming convention changes and
# releases where goreleaser may not have run.

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
  echo "Known working versions: v1.30.0, v1.28.3" >&2
  exit 1
fi

# Find the asset matching our OS and arch (case-insensitive for robustness)
DOWNLOAD_URL=""
while IFS= read -r url; do
  if echo "$url" | grep -qi "${OS}" && echo "$url" | grep -qi "${ARCH}"; then
    DOWNLOAD_URL="$url"
    break
  fi
done <<< "$ASSET_URLS"

if [[ -z "$DOWNLOAD_URL" ]]; then
  echo "Error: no asset found for ${OS}/${ARCH} in release ${TAG_NAME}." >&2
  echo "" >&2
  echo "Available assets:" >&2
  echo "$ASSET_URLS" | sed 's/^/  /' >&2
  exit 1
fi

ARCHIVE_NAME="$(basename "$DOWNLOAD_URL")"
echo "Asset:        ${ARCHIVE_NAME}"
echo "Download URL: ${DOWNLOAD_URL}"

# ── download ────────────────────────────────────────────────────────────────
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

ARCHIVE_PATH="${TMPDIR}/${ARCHIVE_NAME}"

echo "Downloading..."
http_download "$DOWNLOAD_URL" "$ARCHIVE_PATH"

# ── extract ─────────────────────────────────────────────────────────────────
echo "Extracting..."
if [[ "$ARCHIVE_NAME" == *.zip ]]; then
  if command -v unzip &>/dev/null; then
    unzip -q -o "$ARCHIVE_PATH" -d "$TMPDIR"
  elif command -v powershell &>/dev/null; then
    powershell -Command "Expand-Archive -Path '$ARCHIVE_PATH' -DestinationPath '$TMPDIR' -Force"
  else
    echo "Error: unzip or powershell is required to extract .zip files." >&2
    exit 1
  fi
elif [[ "$ARCHIVE_NAME" == *.tar.gz || "$ARCHIVE_NAME" == *.tgz ]]; then
  tar -xzf "$ARCHIVE_PATH" -C "$TMPDIR"
else
  echo "Error: unknown archive format: ${ARCHIVE_NAME}" >&2
  exit 1
fi

# ── locate the binary ──────────────────────────────────────────────────────
if [[ "$OS" == "Windows" ]]; then
  BINARY_NAME="pmtiles.exe"
else
  BINARY_NAME="pmtiles"
fi

EXTRACTED_BINARY="${TMPDIR}/${BINARY_NAME}"

# Some archives may nest the binary in a subdirectory — search for it
if [[ ! -f "$EXTRACTED_BINARY" ]]; then
  EXTRACTED_BINARY="$(find "$TMPDIR" -name "$BINARY_NAME" -type f | head -1)"
fi

if [[ -z "$EXTRACTED_BINARY" || ! -f "$EXTRACTED_BINARY" ]]; then
  echo "Error: could not find '${BINARY_NAME}' in the extracted archive." >&2
  echo "Archive contents:" >&2
  find "$TMPDIR" -type f | sed 's/^/  /' >&2
  exit 1
fi

# ── install ─────────────────────────────────────────────────────────────────
mkdir -p "$DEST"
mv "$EXTRACTED_BINARY" "${DEST}/${BINARY_NAME}"
chmod +x "${DEST}/${BINARY_NAME}"

INSTALLED_PATH="$(cd "$DEST" && pwd)/${BINARY_NAME}"
echo ""
echo "Installed: ${INSTALLED_PATH}"
echo "Version:   $("${INSTALLED_PATH}" version 2>/dev/null || echo "${TAG_NAME}")"
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
