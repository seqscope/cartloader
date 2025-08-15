import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict
import tarfile
from cartloader.utils.utils import write_dict_to_file, load_file_to_dict
from cartloader.utils.file_helper import find_valid_path, find_valid_path_from_zip

## =====INFO: IMGs===
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


def detect_xenium_output(_args):
    # Categorized argument groups
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Detect the 10X Xenium Ranger output files and summarize their locations into a json file.")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
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
        datdict=load_file_to_dict(args.out_json)
    else:
        datdict={}
    
    for key, pattern in xenium_key2patterns.items():
        if datdict.get(key, None) is None:
            datdict[key] = find_valid_path(pattern, args.in_dir)
            if datdict[key] is None and len(pattern.get("zips", []))>0:
                print("test")
                datdict[key] = find_valid_path_from_zip(pattern, args.in_dir, args.unzip_dir, args.overwrite)
        # Validate required file
        if pattern.get("required", False) and datdict[key] is None:
            filenames = ",".join(pattern.get("filename", []))
            raise ValueError(f"Cannot find the required file for '{key}' using filename pattern(s): {filenames}")
    
    # drop the pairs with None values
    datdict={k: v for k, v in datdict.items() if v is not None}
    write_dict_to_file(datdict, args.out_json, check_equal=True, sort_keys=False)
 
if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
