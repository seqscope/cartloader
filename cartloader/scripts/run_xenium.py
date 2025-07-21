import sys, os, argparse, logging, subprocess
import pandas as pd
from pathlib import Path
from cartloader.utils.minimake import minimake

from cartloader.utils.utils import add_param_to_cmd, cmd_separator

def parse_arguments(_args):
    """
    Build resources for CartoScope, joining the pixel-level results from FICTURE
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader run_cartload2", description="Build resources for CartoScope, joining the pixel-level results from FICTURE")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--load-xeniumranger', action='store_true', default=False, help='Load the xenium data from a directory, and summarize paths of xenium files into a json file')
    cmd_params.add_argument('--import-cells', action='store_true', default=False, help='')
    cmd_params.add_argument('--run-fig2pmtiles', action='store_true', default=False, help='')
    cmd_params.add_argument('--run-cartload2', action='store_true', default=False, help='')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/output directory/files. --cell-dir is optional. --fic-dir ")
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output directory')
    inout_params.add_argument('--xenium-dir', type=str, help='Path to input Xenium directory containing its cell & boundary & histology output files')
    inout_params.add_argument('--sge-dir', type=str, help='Path to input SGE directory from sge_convert. It should contains SGE and a YAML/JSON file summarizing the path of transcript, feature, and minmax. The file name of this YAML/JSON file is provided via --in-sge-assets. If --fic-dir is provided, --sge-dir can be skipped')
    inout_params.add_argument('--fic-dir', type=str, help='Path to input FICTURE directory containing FICTURE output. It should contains a YAML/JSON file file provding the SGE and the FICTURE parameters.The file name of this YAML/JSON file is specified by --in-fic-params.')
    inout_params.add_argument('--out-cell', type=str, help='Output stem (i.e., prefix without path) of cell and boundary')
    inout_params.add_argument('--ome-tif', type=str, help="Path to the input ome tiff.")

    key_params = parser.add_argument_group("Key Parameters", "Key parameters frequently used by users")
    key_params.add_argument('--id', type=str, required=True, help='The identifier of the output assets')
    key_params.add_argument('--title', type=str, help='The title of the output assets')
    key_params.add_argument('--desc', type=str, help='The description of output assets')
    key_params.add_argument('--cell-id', type=str, help='Identifier of the cell factor (defaults to --out-cell)')
    key_params.add_argument('--cell-name', type=str, help='Name of the cell factor (defaults to --out-cell)')
    # aux_hist_params = parser.add_argument_group("Auxiliary fig2pmtiles Parameters (using default is recommended)")

    aux_xenium_params = parser.add_argument_group("Auxiliary Import Xenium Parameters. If --load-xenium is enabled, cartloader will detect those files ")
    aux_xenium_params.add_argument('--csv-cells', type=str, help='Name of the CSV file containing the cell locations')
    aux_xenium_params.add_argument('--csv-boundaries', type=str, help='Name of the CSV file containing the cell boundary files')
    aux_xenium_params.add_argument('--csv-clust', type=str, help='Name of the CSV file containing the cell clusters')
    aux_xenium_params.add_argument('--csv-diffexp', type=str, help='Name of the CSV file containing the differential expression results')
    aux_xenium_params.add_argument('--tsv-cmap', type=str, help='Name of the TSV file containing the color map for the cell clusters')
    aux_xenium_params.add_argument('--de-max-pval', type=float, default=0.01, help='Maximum p-value threshold for differential expression')
    aux_xenium_params.add_argument('--de-min-fc', type=float, default=1.2, help='Minimum fold change threshold for differential expression')

    aux_cartload_params = parser.add_argument_group("Auxiliary run_cartload2 Parameters (using default is recommended)")
    aux_cartload_params.add_argument('--in-sge-assets', type=str, default="sge_assets.json", help='The YAML/JSON file containing both the SGE files')
    aux_cartload_params.add_argument('--in-fic-params', type=str, default="ficture.params.json", help='The YAML/JSON file containing both the SGE files and FICTURE parameters')
    aux_cartload_params.add_argument('--colname-feature', type=str, default='gene', help='Input/output Column name for gene name (default: gene)')
    aux_cartload_params.add_argument('--colname-count', type=str, default='count', help='Column name for feature counts')

    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gzip', type=str, default="gzip", help='Path to gzip binary. For faster processing, use "pigz -p4"')
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    env_params.add_argument('--tippecanoe', type=str, default=f"tippecanoe", help='Path to tippecanoe binary') # default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", 
    env_params.add_argument('--spatula', type=str, default=f"spatula",  help='Path to spatula binary') # default=f"{repo_dir}/submodules/spatula/bin/spatula",

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)

xenium_argset=[
    "csv_cells", "csv_boundaries", "csv_clust", "csv_diffexp", "tsv_cmap", 
    "de_max_pval", "de_min_fc",
    "tippecanoe"
]

cartload_argset =[
    "id", "title", "desc",
    "in_sge_assets", "in_fic_params","colname_feature","colname_count",
    "gzip", "pmtiles", "gdaladdo", "tippecanoe", "spatula"
]

def run_xenium(_args):
    """
    Build resources for CartoScope
    """
    # parse argument
    args=parse_arguments(_args)
    
    # mm
    mm = minimake()

    xenium_file_json=f"{args.out_dir}/{args.out_cell}.json"
    if args.load_xenium:
        if args.cell_auto:
            # summarize the 
            cmds = cmd_separator([], f"Create a JSON file summarizing Xenium output from {args.xenium_dir}")
            load_xenium_cmd=" ".join([
                "cartloader", "detect_xenium_output",
                f"--in-dir {args.xenium_dir}"
                f"--out-json {xenium_file_json}",
                f"--unzip-dir {args.out_dir}/unzip",
            ])
            cmds.append(load_xenium_cmd)
            mm.add_target(xenium_file_json, [], cmds) 

    if args.import_cells:
        import_cell_cmd=" ".join([
            "cartloader", "import_xenium_output",
            "--all",
            f"--indir {args.cell_dir}"
            f"--outprefix {args.out_dir}/{args.out_cell}",
            f"--id {args.cell_id}" if args.cell_id else "",
            f"--name {args.cell_name}" if args.cell_name else "",
        ])
    
        import_cell_cmd = add_param_to_cmd(import_cell_cmd, args, xenium_argset)
        import_cell_cmd = add_param_to_cmd(import_cell_cmd, args, ["restart","n_jobs","threads"])
        cmds.append(import_cell_cmd)

    if args.run_cartload2:
        cartload_cmd = " ".join([
            "cartloader", "run_cartload2_generic",
            f"--out-dir {args.out_dir}",
            f"--sge_dir {args.sge_dir}" if args.sge_dir else "",
            f"--fic-dir {args.fic_dir}" if args.fic_dir else "",
            f"--cell-assets {args.out_dir}/{args.out_cell}" if args.import_xenium else "",
            # f"--background-assest {args.out_dir}/{args.out_cell}"
        ])
        cartload_cmd = add_param_to_cmd(cartload_cmd, args, cartload_argset)
        cartload_cmd = add_param_to_cmd(cartload_cmd, args, ["restart","n_jobs","threads"])
        cmds.append(cartload_cmd)

    for cmd in cmds:
        print(cmd)
        if not args.dry_run:
            result = subprocess.run(cmd, shell=True)
            if result.returncode != 0:
                print(f"Error in executing: {cmd}")
                sys.exit(1)


# cartloader run_fig2pmtiles \
#     --makefn run_fig2pmtiles_morphology-focus.ome.mk \
#     --transform \
#     --upper-thres-quantile 0.95 \
#     --lower-thres-quantile 0.5   \
#     --colorize 0000FF \
#     --georeference \
#     --geotif2mbtiles \
#     --mbtiles2pmtiles \
#     --update-catalog \
#     --basemap-key HnE_focus \
#     --in-fig /net/1000g/hmkang/data/xenium/ff_mouse_brain/section1/morphology_focus.ome.tif \
#     --out-prefix /net/1000g/hmkang/weiqiuc/cart/xenium-mouse-brain-ff-section1/xenium-mouse-brain-ff-section1-v1/cartload/morphology-focus.ome      \
#     --n-jobs 10  \
#     --pmtiles /net/1000g/hmkang/weiqiuc/tools/bin/pmtiles

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
