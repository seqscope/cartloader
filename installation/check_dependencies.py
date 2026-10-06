
import subprocess
import sys
import shutil
import importlib.util
import importlib.metadata
import os
import re
from pathlib import Path

# Identify Repo Root
# script is in <root>/installation/check_dependencies.py
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

def check_command(cmd, name=None, fallback_paths=None):
    if name is None:
        name = cmd
    
    # 1. Check PATH
    path = shutil.which(cmd)
    if path:
        print(f"[OK] {name:<20} found at {path}")
        return True
    
    # 2. Check Fallbacks
    if fallback_paths:
        for p in fallback_paths:
            full_path = REPO_ROOT / p
            if full_path.exists() and os.access(full_path, os.X_OK):
                print(f"[OK] {name:<20} found at {full_path} (local build)")
                return True
            # Also check if it's just a directory (e.g. for some tools) or if the binary is inside
            
    print(f"[MISSING] {name:<20} not found in PATH or local build locations")
    return False

def check_submodule(submodule_path):
    full_path = REPO_ROOT / submodule_path
    if full_path.exists() and any(full_path.iterdir()):
        print(f"[OK] Submodule {submodule_path} appears initialized")
        return True
    else:
        print(f"[MISSING] Submodule {submodule_path} is missing or empty")
        return False

# pip distribution name -> import name, where they differ beyond '-' -> '_'
IMPORT_NAMES = {
    "Pillow": "PIL",
    "PyYAML": "yaml",
    "scikit-learn": "sklearn",
    "google-genai": "google.genai",
    "opentsne": "openTSNE",
}

def read_requirements():
    """(name, minimum version or None, extra or None) for every dependency the package
    declares, read from pyproject.toml so this check always matches what the install step
    installs. Python < 3.11 has no tomllib; there the installed package's metadata is read
    instead (importlib.metadata.PackageNotFoundError if cartloader is not installed)."""
    try:
        import tomllib
        with open(REPO_ROOT / "pyproject.toml", "rb") as f:
            project = tomllib.load(f)["project"]
        specs = [(s, None) for s in project.get("dependencies", [])]
        for extra, deps in project.get("optional-dependencies", {}).items():
            specs += [(s, extra) for s in deps]
    except ModuleNotFoundError:
        specs = []
        for s in importlib.metadata.requires("cartloader") or []:
            m = re.search(r"""extra\s*==\s*['"]([^'"]+)['"]""", s)
            specs.append((s.split(";", 1)[0], m.group(1) if m else None))
    reqs = []
    for spec, extra in specs:
        name = re.match(r"[A-Za-z0-9][A-Za-z0-9._-]*", spec.strip()).group(0)
        m = re.search(r">=\s*([0-9][0-9.]*)", spec)
        reqs.append((name, m.group(1) if m else None, extra))
    return reqs

def _version_tuple(v):
    return tuple(int(p) for p in re.findall(r"\d+", v)[:3])

def check_python_module(module_name, min_version=None, extra=None):
    if module_name == "parquet-tools":
        # Skip package check, handled as binary
        return True
    import_name = IMPORT_NAMES.get(module_name, module_name.replace("-", "_"))
    note = f" (optional: [{extra}])" if extra else ""

    try:
        if importlib.util.find_spec(import_name) is not None:
            if min_version:
                installed = importlib.metadata.version(module_name)
                if _version_tuple(installed) < _version_tuple(min_version):
                    print(f"[OUTDATED] {module_name:<20} {installed} installed, >= {min_version} required{note}")
                    return False
            print(f"[OK] {module_name:<20} installed{note}")
            return True
    except Exception:
        pass

    print(f"[MISSING] {module_name:<20} not installed{note}")
    return False

def check_r_package(package_name):
    cmd = ["Rscript", "-e", f'if (!requireNamespace("{package_name}", quietly = TRUE)) quit(status = 1)']
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"[OK] {package_name:<20} installed")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"[MISSING] {package_name:<20} not installed")
        return False

def main():
    print(f"Repo Root detected as: {REPO_ROOT}\n")

    # Check for globally available tools first to inform submodule checks
    has_imagemagick = shutil.which("convert") is not None or shutil.which("magick") is not None

    print("Checking Submodules...")
    # List of submodules to check
    submodules = [
        "submodules/spatula",
        "submodules/tippecanoe",
        "submodules/punkst",
        "submodules/ImageMagick", # Conditional
        "submodules/ficture"
    ]
    missing_submodules = []
    
    for sub in submodules:
        if sub == "submodules/ImageMagick" and has_imagemagick:
            print(f"[NOTE] Skipping {sub} check because ImageMagick is in PATH")
            continue
            
        if not check_submodule(sub):
            missing_submodules.append(sub)

    print("\nChecking System Utilities & Binaries...")
    # Map tool name to potential fallback local paths (relative to REPO_ROOT)
    # Binaries that we might build locally
    tool_fallbacks = {
        "spatula": ["submodules/spatula/build/spatula", "submodules/spatula/bin/spatula"],
        "tippecanoe": ["submodules/tippecanoe/tippecanoe"],
        "pmtiles": [], # User installed manually
    }
    
    system_tools = [
        "gzip", "sort", "bc", "perl", 
        "bgzip", "tabix",             
        "aws",                        
        "gdal_translate", "gdaladdo", "gdalinfo", 
        "spatula",                    
        "tippecanoe",                 
        "pmtiles",                    
        "Rscript",                    
        "python"                      
    ]
    
    missing_tools = []
    for tool in system_tools:
        fallbacks = tool_fallbacks.get(tool)
        if not check_command(tool, fallback_paths=fallbacks):
            missing_tools.append(tool)

    # ImageMagick (used by spatula's CImg): ImageMagick 6 provides `convert` (e.g. Ubuntu's
    # apt package) and ImageMagick 7 provides `magick`; either one will do.
    if shutil.which("convert"):
        check_command("convert", name="ImageMagick")
    elif not check_command("magick", name="ImageMagick",
                           # depends on how IM is built, usually make install puts it in /usr/local
                           fallback_paths=["submodules/ImageMagick/utilities/magick"]):
        missing_tools.append("ImageMagick (convert or magick)")
            
    print("\nChecking Python Packages...")
    missing_python, missing_optional = [], {}
    try:
        requirements = read_requirements()
    except importlib.metadata.PackageNotFoundError:
        print("[MISSING] cartloader is not installed, so its dependencies cannot be read "
              "(install it with: pip install -e .)")
        requirements = []
        missing_python.append("cartloader")
    for pkg, min_version, extra in requirements:
        if not check_python_module(pkg, min_version, extra):
            name = pkg if min_version is None else f"{pkg}>={min_version}"
            if extra is None:
                missing_python.append(name)
            else:
                missing_optional.setdefault(extra, []).append(name)

    # Check parquet-tools CLI specifically
    if not check_command("parquet-tools"):
         missing_tools.append("parquet-tools")

    print("\nChecking R Packages...")
    r_packages = ["argparse", "data.table", "RcppParallel", "ggplot2", "uwot"]
    
    missing_r = []
    if check_command("Rscript"):
        for pkg in r_packages:
            if not check_r_package(pkg):
                missing_r.append(pkg)
    else:
        print("Skipping R package check because Rscript is missing.")
        missing_r = r_packages

    print("\n" + "="*40)
    print("Summary")
    print("="*40)
    
    clean = True
    if missing_submodules:
        print(f"Missing/Empty Submodules: {', '.join(missing_submodules)}")
        clean = False
    if missing_tools:
        print(f"Missing System/External Tools: {', '.join(missing_tools)}")
        clean = False
    if missing_python:
        print(f"Missing Python Packages: {', '.join(missing_python)}")
        clean = False
    if missing_r:
        print(f"Missing R Packages: {', '.join(missing_r)}")
        clean = False
    # Optional extras do not fail the check; they are needed only for the features they cover.
    for extra, names in missing_optional.items():
        print(f"Missing optional [{extra}] packages (install with: pip install -e '.[{extra}]'): {', '.join(names)}")

    if clean:
        print("All dependencies appear to be properly installed!")
        sys.exit(0)
    else:
        print("Some dependencies are missing. Please verify your installation.")
        sys.exit(1)

if __name__ == "__main__":
    main()
