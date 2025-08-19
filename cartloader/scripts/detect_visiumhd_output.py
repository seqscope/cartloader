import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict
import tarfile
from pathlib import Path

from cartloader.utils.utils import write_dict_to_file, load_file_to_dict
from cartloader.utils.file_helper import find_valid_path, find_valid_path_from_zip

## =====INFO=========
# (source: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview)
# barcode_mappings.parquet	This file efficiently stores spatial mapping information, essentially functioning as a CSV that tracks relationships between barcodes (squares), nuclei, cells, and bins within your Visium HD data. See the Segmented Outputs page for more details.
# binned_outputs	        By default, this directory has three subdirectories: square_002um, square_008um, and square_016um. Each directory contains filtered_feature_bc_matrix, raw_feature_bc_matrix, spatial, filtered_feature_bc_matrix.h5, and raw_feature_bc_matrix.h5. The analysis directory is only provided at 8 and 16 µm bin size. The cloupe.cloupe is only provided at 8 µm bin size. The raw_probe_bc_matrix.h5 is only provided at 2 µm resolution.
# cloupe_008um.cloupe	    Symlink to the .cloupe file at 8 µm bin size
# cloupe_cell.cloupe	    Symlink to the .cloupe file with cell segmentation
# feature_slice.h5	        A new file type, specific to Visium HD, to support efficient fetching of 2 µm resolution image slices for a single gene or multiple genes. See this page for more details.
# metrics_summary.csv	    Run summary metrics in CSV format
# molecule_info.h5	        Contains per-molecule information for all molecules that contain a valid barcode, valid UMI, and were assigned with high confidence to a gene barcode or bin.
# probe_set.csv	            Copy of the input probe set reference CSV file.
# segmented_outputs	        Folder containing segmented outputs. Contains analysis, cell_segmentations.geojson, cloupe.cloupe, filtered_feature_cell_matrix, filtered_feature_cell_matrix.h5, graphclust_annotated_cell_segmentations.geojson, graphclust_annotated_nucleus_segmentations.geojson, nucleus_segmentations.geojson, raw_feature_cell_matrix, raw_feature_cell_matrix.h5, and spatial. See the Segmented Outputs page for more details.
# spatial	                Folder containing outputs that capture the spatiality of the data. See the Spatial Outputs page for more details.
# web_summary.html	        Run summary metrics and plots in HTML format
## ==============

def populate_from_suffixes(pattern: dict, in_dir: str, dest_key: str, suffix_key: str) -> None:
    """
    If pattern[dest_key] is missing and pattern[suffix_key] exists,
    glob for *<suffix> files in in_dir and set pattern[dest_key] to a
    list of basenames (order-preserving, de-duplicated).
    """
    # Already populated or no suffixes to search for
    if pattern.get(dest_key) or not pattern.get(suffix_key):
        return

    suffixes = pattern[suffix_key]
    if isinstance(suffixes, str):
        suffixes = [suffixes]

    matches = []
    in_path = Path(in_dir)
    for sfx in suffixes:
        for p in in_path.glob(f"*{sfx}"):
            if p.is_file():
                matches.append(p.name)

    # De-duplicate while preserving order
    deduped = list(dict.fromkeys(matches))
    pattern[dest_key] = deduped

# ? square_008um and 16 contains diffexp and clustering, pca, umap
visiumhd_key2patterns={
    # SGE: 
    # use filtered_feature_bc_matrix, which contains only tissue-associated barcodes,
    "TRANSCRIPT_MEX": {
        "required": True,
        "filenames":["square_002um/filtered_feature_bc_matrix", "binned_outputs/square_002um/filtered_feature_bc_matrix"],
        "zip_suffixes": ["_square_002um_binned_outputs.tar.gz", "_binned_outputs.tar.gz"]
    },
    "SCALE": {
        "required": True,
        "filenames":["square_002um/spatial/scalefactors_json.json", "binned_outputs/square_002um/spatial/scalefactors_json.json"],
        "zip_suffixes": ["_square_002um_binned_outputs.tar.gz", "_binned_outputs.tar.gz"]
    },
    "POSITION":{
        "required": True,
        "filenames":["square_002um/spatial/tissue_positions.parquet", "binned_outputs/square_002um/spatial/tissue_positions.parquet"],
        "zip_suffixes": ["_square_002um_binned_outputs.tar.gz", "_binned_outputs.tar.gz"]
    },
    # CELL: 
    "CELL_FEATURE_MEX": {
        "required": False,
        "filenames":["segmented_outputs/filtered_feature_cell_matrix"],
        "zip_suffixes": ["_segmented_outputs.tar.gz"]
        },
    "CELL_GEOJSON":{
        "required": False,
        "filenames":["segmented_outputs/cell_segmentations.geojson"],
        "zip_suffixes": ["_segmented_outputs.tar.gz"]
    },
    # "CELL_GEOJSON_ANNOTATED":{
    #     "required": False,
    #     "filenames":["segmented_outputs/graphclust_annotated_cell_segmentations.geojson"],
    #     "zip_suffixes": ["_segmented_outputs.tar.gz"]
    # },
    # * Cluster
    # Barcode,Cluster
    # cellid_000000001-1,22
    "CLUSTER": {
        "required": False,
        "filenames":["segmented_outputs/analysis/clustering/gene_expression_graphclust/clusters.csv"],
        "zip_suffixes": ["_segmented_outputs.tar.gz"]
        },
    "DE":{
        "required": False,
        "filenames":["segmented_outputs/analysis/diffexp/gene_expression_graphclust/differential_expression.csv"],
        "zip_suffixes": ["_segmented_outputs.tar.gz"]            
        },
    # 
    # IMGs: NOT required
    # "HIGH_RESOLUTION_PNG": {
    #     "required": False,
    #     "filenames":["square_002um/filtered_feature_bc_matrix/tissue_hires_image.png", "binned_outputs/square_002um/filtered_feature_bc_matrix/tissue_hires_image.png"],
    #     "zip_suffixes": ["_square_002um_binned_outputs.tar.gz", "_binned_outputs.tar.gz"]
    #     }
    "HnE_BTF": {
        "required": False,
        "filename_suffixes":["_tissue_image.btf"],
        "zip_suffixes": ["_square_002um_binned_outputs.tar.gz", "_binned_outputs.tar.gz"]
        }

}

def detect_visiumhd_output(_args):
    # Categorized argument groups
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Detect the 10X Visium HD output files and summarize their locations into a json file.")
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files.")
    inout_params.add_argument('--in-dir', type=str, required=True, help='Path to input directory where Visium HD files are located.')
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

    # Update the pattern if "zips" is empty but "zip_suffixes" are defined; as well as filenames
    for key, pattern in visiumhd_key2patterns.items():
        populate_from_suffixes(pattern, args.in_dir, dest_key="zips",     suffix_key="zip_suffixes")
        populate_from_suffixes(pattern, args.in_dir, dest_key="filenames", suffix_key="filename_suffixes")

    for key, pattern in visiumhd_key2patterns.items():
        if datdict.get(key, None) is None:
            datdict[key] = find_valid_path(pattern, args.in_dir)
            if datdict[key] is None and len(pattern.get("zips", []))>0:
                datdict[key] = find_valid_path_from_zip(pattern, args.in_dir, args.unzip_dir, args.overwrite)
                
        # Validate required file
        if pattern.get("required", False) and datdict[key] is None:
            filenames = ",".join(pattern.get("filenames", []))
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
