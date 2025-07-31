import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict
import tarfile
from cartloader.utils.utils import write_dict_to_file, load_file_to_dict


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

## =====INFO===
# Morphology images
# A series of tissue morphology images are output by the pipeline, which are either nuclei-stained (DAPI) or nuclei and multi-tissue stained (DAPI, cell boundary, interior stains) images in OME-TIFF format. These files include a pyramid of resolutions and tiled chunks of image data, which allows for efficient interactive image visualization (JPEG-2000 compression, 16-bit grayscale, full and downsampled resolutions down to 256 x 256 pixels, learn more here). All morphology image files can be read by Xenium Explorer.


# The morphology.ome.tif is a 3D Z-stack of the DAPI image that can be useful to resegment cells, assess segmentation quality, and view data. DAPI image processing is described here.
# The morphology_focus/ directory contains the 2D autofocus projection images for the nuclei DAPI stain image, as well as three additional stain images for Xenium outputs generated with the multimodal cell segmentation assay workflow. These files are in multi-file OME-TIFF format. They each contain a pyramid of images including full resolution and downsampled images. The image order is specified in OME-XML metadata.
# morphology_focus_0000.ome.tif: DAPI image
# morphology_focus_0001.ome.tif: boundary (ATP1A1/E-Cadherin/CD45) image
# morphology_focus_0002.ome.tif: interior - RNA (18S) image
# morphology_focus_0003.ome.tif: interior - protein (alphaSMA/Vimentin) image

# morphology_mip.ome.tif: Cell morphology DAPI with the maximum intensity projection (MIP) of the DAPI Z-stack image.

# he_image.ome.tiff: (NOT FOUND in available datasets) An H&E image of the tissue section, taken after the Xenium protocol, if available. Note, this image is acquired on a different microscope, and is not registered to the morphology image.
## =====

xenium_key2patterns={
    # SGE:
    "TRANSCRIPT": {
        "required": True,
        "filename":["transcripts.csv.gz", "transcripts.parquet"]
    },
    # CELL:
    "CELL": {
        "required": True,
            "filename": ["cells.csv.gz", "cells.parquet"],
        },
    "BOUNDARY":  {
        "required": True,
        "filename": ["cell_boundaries.csv.gz"],
        },
    "CELL_FEATURE_MEX": {
        "required": True,
        "filename":["cell_feature_matrix"],
        "zips": ["cell_feature_matrix.tar.gz"]
        },
    # CLUSTER and DE
    "CLUSTER": {
        "required": True,
        "filename":["analysis/clustering/gene_expression_graphclust/clusters.csv"],
        "zips": ["analysis.tar.gz"]
        },
    "DE":{
        "required": True,
        "filename":["analysis/diffexp/gene_expression_graphclust/differential_expression.csv"],
        "zips": ["analysis.tar.gz"]
        },
    # IMGs: NOT required
    "DAPI_OME": {
        "required": False,
        "filename":["morphology_focus.ome.tif", "morphology_focus/morphology_focus_0000.ome.tif" ]
        },
    "BOUNDARY_OME": {
        "required": False,
        "filename":["morphology_focus/morphology_focus_0001.ome.tif" ]
        },
    "INTERIOR_RNA_OME": {
        "required": False,
        "filename":["morphology_focus/morphology_focus_0002.ome.tif" ]
        },
    "INTERIOR_PROTEIN_OME": {
        "required": False,
        "filename":["morphology_focus/morphology_focus_0003.ome.tif" ]
        },
    "DAPI_3D_OME": {
        "required": False,
        "filename":["morphology.ome.tif"]
        },
    "DAPI_MIP_OME": {
        "required": False,
        "filename":["morphology_mip.ome.tif"]
    }
        
}

def find_valid_path(pattern, key, in_dir, cellbounds):
    for filename in pattern["filename"]:
        f=os.path.join(in_dir,filename)
        if os.path.exists(f):
            cellbounds[key]=f
    return cellbounds

def find_valid_path_from_zip(pattern, key, in_dir, unzip_dir, cellbounds, overwrite=False):
    for zip_fn in pattern["zips"]:
        zip_in=os.path.join(in_dir, zip_fn)
        assert os.path.exists(zip_in), f"An input compressed file does not exist: {zip_in}"

        zip_stem, zip_suffix = split_prefix_suffix_from_compression(zip_fn)

        # decompress the file (unzip_dir the root directory to host the decompressed files)
        unzip_subdir=os.path.join(unzip_dir, zip_stem)
        if not os.path.exists(unzip_subdir) or overwrite:
            if zip_suffix == ".tar.gz":
                decompress_tar_gz(zip_in, unzip_dir)
            else:
                raise ValueError(f"Found a compressed input file with the {zip_suffix} extension, which is not supported by cartloader. Please unzip the file manually before proceeding: {zip_in}")
        if os.path.exists(unzip_subdir) and overwrite:
            print(f"Warning: the decompressed files will be overwritten: {unzip_subdir}")    
        
        cellbounds=find_valid_path(pattern, key, unzip_dir, cellbounds)

    return cellbounds

def detect_xenium_output(_args):
    # Categorized argument groups
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Detect the 10X Xenium Ranger output files and summarize their locations into a json file.")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files. --cell-dir is optional. --fic-dir ")
    inout_params.add_argument('--in-dir', type=str, required=True, help='Path to input directory where Xenium files are located.')
    inout_params.add_argument('--out-json', type=str, required=True, help='Path to the output JSON file that will list all detected files.')
    inout_params.add_argument('--unzip-dir', type=str, help='Optional directory where compressed files from --in-dir will be decompressed (defaults to --in-dir if not specified; ensure write access)')
    inout_params.add_argument('--overwrite',  action='store_true', default=False, help=' If set, existing JSON and decompressed files will be overwritten; otherwise, existing entries and decompressed directories will be reused or skipped.')
    args = parser.parse_args(_args)

    assert os.path.isdir(args.in_dir), f"The --in-dir should be a directory: {args.in_dir}"

    out_dir=os.path.dirname(args.out_json)
    os.makedirs(out_dir, exist_ok=True)

    if args.unzip_dir is None:
        args.unzip_dir=args.in_dir
    os.makedirs(args.unzip_dir, exist_ok=True)

    if os.path.exists(args.out_json) and not args.overwrite:
        cellbounds=load_file_to_dict(args.out_json)
    else:
        cellbounds={}
    
    for key, pattern in xenium_key2patterns.items():
        if cellbounds.get(key, None) is None:
            if pattern["required"]:
                cellbounds[key] = None
            cellbounds = find_valid_path(pattern, key, args.in_dir, cellbounds)
            if pattern["required"] and cellbounds[key] is None and len(pattern.get("zips", []))>0:
                cellbounds = find_valid_path_from_zip(pattern, key, args.in_dir, args.unzip_dir, cellbounds, args.overwrite)

    # arrange the order 
    if all_items_defined(cellbounds):
        write_dict_to_file(cellbounds, args.out_json, check_equal=True, sort_keys=False)
    else:
        raise ValueError(f"Cannot find all required files from {args.in_dir}\nSee:\n{cellbounds}")
    
if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
