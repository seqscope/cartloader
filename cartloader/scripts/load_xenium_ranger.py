import sys, os, gzip, argparse, logging, warnings, shutil, re, inspect, warnings, glob, json
from collections import defaultdict, OrderedDict
import tarfile
from cartloader.utils.utils import write_dict_to_file, load_file_to_dict
from cartloader.utils.file_helper import resolve_paths_by_pattern

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
    "SGE":{
        "TRANSCRIPT": {
            "required": True,
            "filenames":["transcripts.csv.gz", "transcripts.parquet"]
        },
    },
    # CELL-BASED ANALYSIS:
    "CELLS":{
        # * Cell format:
        # "cell_id","x_centroid","y_centroid","transcript_counts","control_probe_counts","genomic_control_counts","control_codeword_counts","unassigned_codeword_counts","deprecated_codeword_counts","total_counts","cell_area","nucleus_area","nucleus_count","segmentation_method"
        # "aaaagkdm-1",170.85508728027344,2017.2412109375,1,0,0,0,0,0,1,46.285157930105925,NaN,0,"Segmented by boundary stain (ATP1A1+CD45+E-Cadherin)"
        "CELL": {
            "required": True,
                "filenames": ["cells.csv.gz", "cells.parquet"],
            },
        # * Cell boundaries:
        # "cell_id","vertex_x","vertex_y","label_id"
        # "aaaagkdm-1",169.3625,2013.0126,1
        "BOUNDARY":  {
            "required": True,
            "filenames": ["cell_boundaries.csv.gz"],
            },
        "CELL_FEATURE_MEX": {
            "required": False,
            "filenames":["cell_feature_matrix"],
            "zips": ["cell_feature_matrix.tar.gz"]
            },
        # * Cluster:
        # Barcode,Cluster
        # aaaagkdm-1,1
        "CLUSTER": {
            "required": True,
            "filenames":["analysis/clustering/gene_expression_graphclust/clusters.csv"],
            "zips": ["analysis.tar.gz"]
            },
        # * DE:
        # Feature ID,Feature Name,Cluster 1 Mean Counts,Cluster 1 Log2 fold change,Cluster 1 Adjusted p value,Cluster 2 Mean Counts,Cluster 2 Log2 fold change,Cluster 2 Adjusted p value,Cluster 3 Mean Counts,Cluster 3 Log2 fold change,Cluster 3 Adjusted p value,Cluster 4 Mean Counts,Cluster 4 Log2 fold change,Cluster 4 Adjusted p value,Cluster 5 Mean Counts,Cluster 5 Log2 fold change,Cluster 5 Adjusted p value,Cluster 6 Mean Counts,Cluster 6 Log2 fold change,Cluster 6 Adjusted p value,Cluster 7 Mean Counts,Cluster 7 Log2 fold change,Cluster 7 Adjusted p value,Cluster 8 Mean Counts,Cluster 8 Log2 fold change,Cluster 8 Adjusted p value,Cluster 9 Mean Counts,Cluster 9 Log2 fold change,Cluster 9 Adjusted p value,Cluster 10 Mean Counts,Cluster 10 Log2 fold change,Cluster 10 Adjusted p value,Cluster 11 Mean Counts,Cluster 11 Log2 fold change,Cluster 11 Adjusted p value,Cluster 12 Mean Counts,Cluster 12 Log2 fold change,Cluster 12 Adjusted p value,Cluster 13 Mean Counts,Cluster 13 Log2 fold change,Cluster 13 Adjusted p value,Cluster 14 Mean Counts,Cluster 14 Log2 fold change,Cluster 14 Adjusted p value,Cluster 15 Mean Counts,Cluster 15 Log2 fold change,Cluster 15 Adjusted p value,Cluster 16 Mean Counts,Cluster 16 Log2 fold change,Cluster 16 Adjusted p value,Cluster 17 Mean Counts,Cluster 17 Log2 fold change,Cluster 17 Adjusted p value,Cluster 18 Mean Counts,Cluster 18 Log2 fold change,Cluster 18 Adjusted p value,Cluster 19 Mean Counts,Cluster 19 Log2 fold change,Cluster 19 Adjusted p value,Cluster 20 Mean Counts,Cluster 20 Log2 fold change,Cluster 20 Adjusted p value,Cluster 21 Mean Counts,Cluster 21 Log2 fold change,Cluster 21 Adjusted p value
        # ENSG00000166535,A2ML1,0.008082077992052654,2.292429464470433,0.000013831178435443352,0.0015317620832224837,-0.2717658073791789,0.16427980780868393,0.001226843332106492,-0.42095084864144816,0.46975425102829566,0.0037832334823669654,1.1762824711906763,0.00626099448573822,0.0009370316994608665,-0.9931426193918895,0.00012420457106652693,0.00035171591775058056,-1.9689174885349274,0.002636099946030819,0.0018219664759913663,0.0766263739653219,0.8995675333094153,0.0009475336160580633,-0.8343888826514707,0.06727367639920934,0.004361559793501598,1.5130532082572863,0.0000000000007537933872968394,0.000846294895303063,-0.992227102994244,0.028453359285143767,0.0018146574418219585,0.051194655215844875,0.9304123002179088,0.00045951336258476185,-1.7421234527453535,0.0013628094548944181,0.0023563153416359694,0.4774816523131893,0.2965413340700354,0.0021212060341883753,0.30481838952695206,0.4926601159653632,0.0014237107155417482,-0.14271909128282623,0.8500806231615298,0.0012562334587749313,-0.3507612852540696,0.6384951496662045,0.002355140186915884,0.5936824059089432,0.6436592400308225,0,-1.0564604491458525,0.697077988750366,0.0037541549892459126,1.6545227197901715,0.5114958184195594,0.0016831868337385438,0.9092558204592383,0.9557130750504969,0.004746991923520686,1.8259931272624748,0.4816191969175982
        "DE":{
            "required": True,
            "filenames":["analysis/diffexp/gene_expression_graphclust/differential_expression.csv"],
            "zips": ["analysis.tar.gz"]
            },
        "PCA_PROJ":{
            "required": False,
            "filenames":["analysis/pca/gene_expression_10_components/projection.csv"],
            "zips": ["analysis.tar.gz"]
            },
        "PCA_VAR":{
            "required": False,
            "filenames":["analysis/pca/gene_expression_10_components/variance.csv"],
            "zips": ["analysis.tar.gz"]
            },
        "UMAP_PROJ":{
            "required": False,
            "filenames":["analysis/umap/gene_expression_2_components/projection.csv"],
            "zips": ["analysis.tar.gz"]
            },
    },
    # IMGs: NOT required
    "IMAGES":{
        "DAPI_OME": {
            "required": False,
            "filenames":["morphology_focus.ome.tif", "morphology_focus/morphology_focus_0000.ome.tif" ]
            },
        "BOUNDARY_OME": {
            "required": False,
            "filenames":["morphology_focus/morphology_focus_0001.ome.tif" ]
            },
        "INTERIOR_RNA_OME": {
            "required": False,
            "filenames":["morphology_focus/morphology_focus_0002.ome.tif" ]
            },
        "INTERIOR_PROTEIN_OME": {
            "required": False,
            "filenames":["morphology_focus/morphology_focus_0003.ome.tif" ]
            },
        # "DAPI_3D_OME": {
        #     "required": False,
        #     "filenames":["morphology.ome.tif"]
        #     },
        "DAPI_MIP_OME": {
            "required": False,
            "filenames":["morphology_mip.ome.tif"]
        }
    }
}


def load_xenium_ranger(_args):
    # Categorized argument groups
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Detect 10x Xenium Ranger files from Xenium Ranger and write a JSON manifest of file locations."
    )
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output paths.")
    inout_params.add_argument('--in-dir', type=str, required=True, help='Path to input directory containing Xenium Ranger outputs')
    inout_params.add_argument('--out-json', type=str, required=True, help='Path to output JSON manifest to write')
    inout_params.add_argument('--unzip-dir', type=str, help='Path to directory for decompressing archives found under --in-dir (defaults to --in-dir)')
    inout_params.add_argument('--overwrite',  action='store_true', default=False, help='Overwrite existing outputs if present')
    args = parser.parse_args(_args)

    assert os.path.exists(args.in_dir), f"Directory not found: {args.in_dir} (--in-dir)"
    assert os.path.isdir(args.in_dir), f"--in-dir should be a directory: {args.in_dir}"

    out_dir=os.path.dirname(args.out_json)
    os.makedirs(out_dir, exist_ok=True)

    if args.unzip_dir is None:
        args.unzip_dir=args.in_dir
    os.makedirs(args.unzip_dir, exist_ok=True)

    if os.path.exists(args.out_json) and not args.overwrite:
        datdict=load_file_to_dict(args.out_json)
    else:
        datdict={}
    
    datdict = resolve_paths_by_pattern(datdict,
                                        xenium_key2patterns,
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
