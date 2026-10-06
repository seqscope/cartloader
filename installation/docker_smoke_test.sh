#!/bin/bash
# Quick checks of a built cartloader image, run inside the container:
#   docker run --rm --entrypoint bash <image> /app/cartloader/installation/docker_smoke_test.sh
#
# Checks that the dependencies are installed (check_dependencies.py), that the bundled
# binaries resolve all their shared libraries, and that every `cartloader <command>`
# imports and has its entry function. Exits non-zero if anything fails.
set -u
REPO_DIR=${1:-/app/cartloader}
fail=0

echo "== Dependencies (check_dependencies.py)"
python3 "${REPO_DIR}/installation/check_dependencies.py" || fail=1

echo
echo "== Bundled binaries"
for b in punkst spatula pmpoint pmtiles tippecanoe geotiff2pmtiles; do
	path=$(command -v "${b}") || { echo "[MISSING] ${b}"; fail=1; continue; }
	if ldd "${path}" 2>/dev/null | grep -q "not found"; then
		echo "[BROKEN]  ${b}: unresolved shared libraries"
		ldd "${path}" | grep "not found"
		fail=1
	else
		echo "[OK]      ${b} (${path})"
	fi
done

echo
echo "== cartloader commands"
# Mirrors cartloader/cli.py: every scripts/*.py is a command whose module must import and
# define a function of the same name.
python3 - <<'EOF' || fail=1
import importlib, os, sys
import cartloader.scripts as scripts
names = sorted(f[:-3] for f in os.listdir(os.path.dirname(scripts.__file__))
               if f.endswith(".py") and f != "__init__.py")
broken = []
for name in names:
    try:
        module = importlib.import_module(f"cartloader.scripts.{name}")
        if not callable(getattr(module, name, None)):
            broken.append((name, f"module defines no function {name}()"))
    except Exception as e:
        broken.append((name, f"{type(e).__name__}: {e}"))
for name, why in broken:
    print(f"[BROKEN]  cartloader {name}: {why}")
print(f"{len(names) - len(broken)}/{len(names)} commands import cleanly")
sys.exit(1 if broken else 0)
EOF

echo
if [ "${fail}" -eq 0 ]; then echo "SMOKE TEST PASSED"; else echo "SMOKE TEST FAILED"; fi
exit ${fail}
