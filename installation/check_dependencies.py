
import subprocess
import sys
import shutil
import importlib.util
import os
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

def check_python_module(module_name):
    import_name = module_name
    if module_name == "Pillow":
        import_name = "PIL"
    elif module_name == "PyYAML":
        import_name = "yaml"
    elif module_name == "scikit-learn":
        import_name = "sklearn"
    elif module_name == "parquet-tools": 
        # Skip package check, handled as binary
        return True
    
    try:
        if importlib.util.find_spec(import_name) is not None:
             print(f"[OK] {module_name:<20} installed")
             return True
    except Exception:
        pass
    
    print(f"[MISSING] {module_name:<20} not installed")
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
    has_magick = shutil.which("magick") is not None

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
        if sub == "submodules/ImageMagick" and has_magick:
            print(f"[NOTE] Skipping {sub} check because 'magick' is in PATH")
            continue
            
        if not check_submodule(sub):
            missing_submodules.append(sub)

    print("\nChecking System Utilities & Binaries...")
    # Map tool name to potential fallback local paths (relative to REPO_ROOT)
    # Binaries that we might build locally
    tool_fallbacks = {
        "spatula": ["submodules/spatula/build/spatula", "submodules/spatula/bin/spatula"],
        "tippecanoe": ["submodules/tippecanoe/tippecanoe"],
        "magick": ["submodules/ImageMagick/utilities/magick"], # depends on how IM is built, usually make install puts it in /usr/local
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
        "magick",                     
        "Rscript",                    
        "python"                      
    ]
    
    missing_tools = []
    for tool in system_tools:
        fallbacks = tool_fallbacks.get(tool)
        if not check_command(tool, fallback_paths=fallbacks):
            missing_tools.append(tool)
            
    print("\nChecking Python Packages...")
    python_packages = [
        "numpy", "pandas", "ficture", "tifffile", "imagecodecs", 
        "Pillow", "psutil", "rasterio", "requests", "setuptools", 
        "PyYAML", "polars"
    ]
    
    missing_python = []
    for pkg in python_packages:
        if not check_python_module(pkg):
            missing_python.append(pkg)

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
        
    if clean:
        print("All dependencies appear to be properly installed!")
        sys.exit(0)
    else:
        print("Some dependencies are missing. Please verify your installation.")
        sys.exit(1)

if __name__ == "__main__":
    main()
