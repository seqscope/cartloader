import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json

from cartloader.utils.utils import write_dict_to_file, load_file_to_dict
from cartloader.utils.file_helper import resolve_paths_by_pattern, populate_from_suffixes

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


## PCA & UMAP
# │   ├── pca
# │   │   └── gene_expression_10_components
# │   │       ├── components.csv    # PC * gene matrix of weights
# │   │       ├── dispersion.csv    # normalized dispersion (variance / mean) --- used in highly variable gene (HVG) selection 
# │   │       ├── features_selected.csv 
# │   │       ├── projection.csv    # cell/barcode coordinates in PCA space
# │   │       └── variance.csv
# │   └── umap
# │       └── gene_expression_2_components
# │           └── projection.csv

## ==============



# ========== Build visiumhd_key2patterns ========

analysis_map={
    "CLUSTER": "analysis/clustering/gene_expression_graphclust/clusters.csv",
    "DE": "analysis/diffexp/gene_expression_graphclust/differential_expression.csv",
    "PCA_PROJ": "analysis/pca/gene_expression_10_components/projection.csv",
    "PCA_VAR": "analysis/pca/gene_expression_10_components/variance.csv",
    "UMAP_PROJ": "analysis/umap/gene_expression_2_components/projection.csv"
}


cell_map= {
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
    "NUCLEUS_GEOJSON":{
            "required": False,
            "filenames":["segmented_outputs/nucleus_segmentations.geojson"],
            "zip_suffixes": ["_segmented_outputs.tar.gz"]
        },
    }

for key, relpath in analysis_map.items():
    cell_map[key] = {
        "required": True if key in ["CLUSTER", "DE"] else False,
        "filenames": [f"segmented_outputs/{relpath}"],
        "zip_suffixes": ["_segmented_outputs.tar.gz"],
    }


def define_hexmap(r):
    hex_pattern = {
        "HEX_FEATURE_MEX":  "filtered_feature_bc_matrix",
        "SCALE":        "spatial/scalefactors_json.json",
        "POSITION":     "spatial/tissue_positions.parquet",
    }

    hex_pattern.update(analysis_map)

    sq_str = f"square_00{r}um" if r < 10 else f"square_0{r}um"
    hex_map = {}
    for key, val in hex_pattern.items():
        hex_map[key] = {
            "required": False,
            "filenames": [f"{sq_str}/{val}", f"binned_outputs/{sq_str}/{val}"],
            "zip_suffixes": [f"_{sq_str}_binned_outputs.tar.gz", "_binned_outputs.tar.gz"],
        }
    return hex_map


#==================

def load_space_ranger(_args):
    # Categorized argument groups
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Detect 10x Visium HD files from Space Ranger and write a JSON manifest of file locations."
    )
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output paths.")
    inout_params.add_argument('--in-dir', type=str, required=True, help='Path to input directory containing Visium HD Space Ranger outputs')
    inout_params.add_argument('--out-json', type=str, required=True, help='Path to output JSON manifest to write')
    inout_params.add_argument('--unzip-dir', type=str, help='Path to directory for decompressing archives found under --in-dir (defaults to --in-dir)')
    inout_params.add_argument('--overwrite',  action='store_true', default=False, help='Overwrite existing outputs if present (default: False)')
    args = parser.parse_args(_args)

    assert os.path.exists(args.in_dir), f"Directory not found: {args.in_dir} (--in-dir)"
    assert os.path.isdir(args.in_dir), f"--in-dir should be a directory: {args.in_dir}"

    out_dir=os.path.dirname(args.out_json)
    os.makedirs(out_dir, exist_ok=True)

    if args.unzip_dir is None:
        args.unzip_dir=args.in_dir
    os.makedirs(args.unzip_dir, exist_ok=True)

    # Build visiumhd_key2patterns 
    visiumhd_key2patterns={
        # SGE: 
        # use filtered_feature_bc_matrix, which contains only tissue-associated barcodes,
        "SGE":{
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
        },
        "CELLS": cell_map,
        "HEX_8um": define_hexmap(8),
        "HEX_16um": define_hexmap(16),
        "IMAGES":  {
            "HnE_BTF": {
                "required": False,
                "filename_suffixes": ["_tissue_image.btf"],
                "zip_suffixes": ["_square_002um_binned_outputs.tar.gz", "_binned_outputs.tar.gz"],
            }
        }

    }

    # - Update the pattern if "zips" is empty but "zip_suffixes" are defined; as well as filenames
    for group, patterns in visiumhd_key2patterns.items():
        for key, spec in patterns.items():
            populate_from_suffixes(spec, args.in_dir, dest_key="zips",     suffix_key="zip_suffixes")
            populate_from_suffixes(spec, args.in_dir, dest_key="filenames", suffix_key="filename_suffixes")

    # Build datdict, return a nested dict
    if os.path.exists(args.out_json) and not args.overwrite:
        datdict=load_file_to_dict(args.out_json)
    else:
        datdict={}

    datdict = resolve_paths_by_pattern(datdict,
                                        visiumhd_key2patterns,
                                        args.in_dir,
                                        args.unzip_dir,
                                        args.overwrite)
                
    write_dict_to_file(datdict, args.out_json, check_equal=True, sort_keys=False)
    print(f"Created a JSON file at: {args.out_json}")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    print(f"Running {script_name} script")

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
