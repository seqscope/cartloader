import os, tarfile
from pathlib import Path


def split_prefix_suffix_from_compression(filename):
    known_suffixes = [
        '.tar.gz', '.tar.bz2', '.tar.xz', '.tar.zst',
        '.zarr.zip', '.zarr.tar.gz', '.zarr.tar',
        '.zip', '.gz', '.bz2', '.xz', '.zst'
    ]

    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    
    suffix = None
    for known_suffix in sorted(known_suffixes, key=len, reverse=True):
        if basename.endswith(known_suffix):
            suffix = known_suffix
            break

    assert suffix is not None, f"Unknown compression suffix in '{filename}'"

    # remove suffix only from the end of the basename
    prefix_basename = basename[:-len(suffix)]
    prefix = os.path.join(dirname, prefix_basename) if dirname else prefix_basename
    return prefix, suffix 

def decompress_tar_gz(file_path, destination_directory):
    try:
        with tarfile.open(file_path, "r:gz") as tar:
            tar.extractall(destination_directory)
        print(f"Successfully decompressed '{file_path}' to '{destination_directory}'.")
    except tarfile.ReadError:
        print(f"Error: Could not open '{file_path}'. It might not be a valid tar.gz file.")
    except Exception as e:
        print(f"An error occurred during decompression: {e}")

def all_items_defined(d):
    # NOT IN USE 
    """
    Recursively check if all values in a nested dictionary are not None.
    Returns True if all values are defined (no None), otherwise False.
    """
    if isinstance(d, dict):
        for v in d.values():
            if not all_items_defined(v):
                return False
    elif isinstance(d, (list, tuple)):
        for item in d:
            if not all_items_defined(item):
                return False
    else:
        if d is None:
            return False
    return True

def find_valid_path(pattern, in_dir):
    """
    Search for the *first* existing file in 'in_dir' matching the filename pattern
    """
    for filename in pattern["filenames"]:
        f=os.path.join(in_dir, filename)
        if os.path.exists(f):
            return f
    return None

def find_valid_path_from_zip(pattern, in_dir, unzip_dir, overwrite=False):
    """
    Search for the *first* existing file within the compressed file
    """
    assert len(pattern.get("zips", []))>0, "No prebuild zip file name"
    
    for zip_fn in pattern["zips"]:
        zip_in=os.path.join(in_dir, zip_fn)
        assert os.path.exists(zip_in), f"An input compressed file does not exist: {zip_in}"
        zip_stem, zip_suffix = split_prefix_suffix_from_compression(zip_fn)
        unzip_subdir=os.path.join(unzip_dir, zip_stem) # decompress the file (unzip_dir the root directory to host the decompressed files)
        if os.path.exists(unzip_subdir):
            if overwrite:
                print(f"Warning: overwriting decompressed files in: {unzip_subdir}")
        if not os.path.exists(unzip_subdir) or overwrite:
            if zip_suffix == ".tar.gz":
                decompress_tar_gz(zip_in, unzip_dir)
            else:
                raise ValueError(f"Unsupported compressed file extension '{zip_suffix}'. Please unzip manually: {zip_in}")

        f = find_valid_path(pattern, unzip_dir)
        if f is not None:
            return f

    return None


def resolve_paths_by_pattern(
    datdict: dict,
    key2pattern: dict,
    in_dir: str,
    unzip_dir: str,
    overwrite: bool,
) -> dict:
    """
    Resolve file paths for visiumhd_key2patterns into datdict.
    - Fills in missing keys using find_valid_path / find_valid_path_from_zip.
    - Collects missing required files and raises ValueError if any are not found.
    """
    def _resolve_path(spec: dict) -> str | None:
        """Return a valid path for this spec, trying zip sources if needed."""
        path = find_valid_path(spec, in_dir)
        if path is None and spec.get("zips"):
            path = find_valid_path_from_zip(spec, in_dir, unzip_dir, overwrite)
        return path

    errors = []

    for group, patterns in key2pattern.items():
        gdict = datdict.setdefault(group, {})

        for key, spec in patterns.items():
            # Only resolve if missing/None
            if gdict.get(key) is None:
                gdict[key] = _resolve_path(spec)
            
            # If cannot find any path, drop it
            if gdict.get(key) is None:
                gdict.pop(key, None)

            # Collect missing required files for error reporting
            if spec.get("required") and not gdict.get(key):
                filenames = spec.get("filenames") or []
                pretty = ", ".join(filenames) if filenames else "<no patterns provided>"
                errors.append(f"[{group}.{key}] using pattern(s): {pretty}")

    if errors:
        raise ValueError("Cannot find the required file(s):\n  - " + "\n  - ".join(errors))

    return datdict


# === populate the key2pattern from suffixes

def populate_from_suffixes(spec: dict, in_dir: str, dest_key: str, suffix_key: str) -> None:
    """
    If spec[dest_key] is missing and spec[suffix_key] exists,
    glob for *<suffix> files in in_dir and set spec[dest_key] to a
    list of basenames (order-preserving, de-duplicated).
    """
    # Already populated or no suffixes to search for
    if spec.get(dest_key) or not spec.get(suffix_key):
        return

    suffixes = spec[suffix_key]
    if isinstance(suffixes, str):
        suffixes = [suffixes]

    matches: list[str] = []
    in_path = Path(in_dir)
    for sfx in suffixes:
        # search recursively so pre-extracted assets in nested folders are detected
        pattern = f"**/*{sfx}"
        for p in in_path.glob(pattern):
            if p.is_file() or p.is_dir():
                rel_path = p.relative_to(in_path)
                matches.append(str(rel_path))

    # De-duplicate while preserving order
    deduped = list(dict.fromkeys(matches))
    spec[dest_key] = deduped
