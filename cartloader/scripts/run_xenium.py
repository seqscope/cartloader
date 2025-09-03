import sys, os, argparse, logging, subprocess, inspect, re
from pathlib import Path
import hashlib
from typing import Iterable, List, Optional

# from cartloader.utils.minimake import minimake
from cartloader.utils.utils import add_param_to_cmd, cmd_separator, run_command_w_preq, load_file_to_dict, assert_unique
from cartloader.utils.pipeline_helper import validate_general_args, validate_imageid_args, validate_imagecol_args, stage_run_ficture2, stage_run_cartload2, stage_upload_aws, stage_upload_zenodo, stage_import_images


def parse_arguments(_args):
    """
    Process Xenium output
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Process 10X Xenium output from Xenium Ranger, and return rastered and tiled sources for CartoScope (detect → sge_convert → ficture2 → images/cells → cartload2 → uploads")

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run mode. Print all commands without executing them.')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    cmd_params = parser.add_argument_group("Commands")
    cmd_params.add_argument('--load-xenium-ranger', action='store_true', default=False, help='Automatically detect Xenium Ranger output files in --xenium-ranger-dir and save their paths to --xenium-ranger-assets (JSON).')
    cmd_params.add_argument('--sge-convert', action='store_true', default=False, help='Convert SGE files from Xenium format into the cartloader-compatible format.')
    cmd_params.add_argument('--run-ficture2', action='store_true', default=False, help='Run FICTURE analysis on the converted SGE, using the latest punkst version.')
    cmd_params.add_argument('--import-ext-ficture2', action='store_true', default=False, help='Import a set of existing FICTURE results.')
    cmd_params.add_argument('--import-cells', action='store_true', default=False, help='Import cell-based analysis results from Xenium Ranger, and add the asset information to catalog.yaml if exists.')
    cmd_params.add_argument('--import-images', action='store_true', default=False, help='Import images as background, such as DAPI, and add the asset information to catalog.yaml if exists.')
    cmd_params.add_argument('--run-cartload2', action='store_true', default=False, help='Package SGE into raster tiles (and import FICTURE results if --run-ficture2 is enabled). A catalog.yaml file will be created to summarize all assets.')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload output assets into AWS')
    cmd_params.add_argument('--upload-zenodo', action='store_true', default=False, help='Upload output assets into Zenodo')

    inout_params = parser.add_argument_group("Input/Output Parameters", 'cartloader determines input files in two ways:	1. Manual input mode (explicit) — If any argument in "Auxiliary Input/Output File/Directories" are specified with --xenium-ranger-dir; 2. Auto-detection mode (implicit) — If no manual input is given, cartloader reads input file paths from --xenium-ranger-assets, which is created by --load-xenium-ranger.')
    inout_params.add_argument('--xenium-ranger-dir', type=str, help='(Required if manual input mode is enabled or running --load-xenium-ranger) Path to the Xenium Ranger output directory containing transcript, cell, boundary, cluster, and image files.')
    inout_params.add_argument('--out-dir', type=str, help='Path to output root directory. If set, --sge-dir, --fic-dir, --cart-dir, --xenium-ranger-assets default to <out_dir>/sge, <out_dir>/ficture2, <out_dir>/cartload2, and <out_dir>/xenium_ranger_assets.json. To use a custom directory layout, skip --out-dir and define those paths manually.')
    inout_params.add_argument('--ext-fic-dir', type=str, help='(Required if --import-ext-ficture2 is enabled) Path to an external FICTURE directory for loading external FICTURE assets.')

    key_params = parser.add_argument_group("Key Parameters")
    # sge convert
    key_params.add_argument('--units-per-um', type=float, default=None, help='(Parameters for --sge-convert) Coordinate unit per um in the input files (default: 1.00).')  
    key_params.add_argument('--filter-by-density', action='store_true', default=False, help='(Parameters for --sge-convert) Filter SGE from format conversion by density (default: False).')
    key_params.add_argument('--exclude-feature-regex', type=str, default=None, help='(Parameters for --sge-convert) A regex pattern of feature/gene names to be excluded (default: "^(BLANK_|DeprecatedCodeword_|NegCon|UnassignedCodeword_)")')
    # ficture2 parameters
    key_params.add_argument('--width', type=str, default=None, help='(Parameters for --run-ficture2; required) Comma-separated hexagon flat-to-flat widths (in um) for LDA training and projection (default: None)') ## Same width will be used for both train width and projection width
    key_params.add_argument('--n-factor', type=str, default=None, help='(Parameters for --run-ficture2; required) Comma-separated list of factor counts for LDA training.')
    key_params.add_argument('--colname-feature', type=str, default='gene', help='(Parameters for --run-ficture2 or --run-cartload2) Input/output column name for gene name (default: gene)')
    key_params.add_argument('--colname-count', type=str, default='count', help='(Parameters for --run-ficture2 or --run-cartload2) Column name for feature counts (default: count)')
    # run_cartload2 parameters
    key_params.add_argument('--id', type=str, help='(Parameters for --run-cartload2; required) Identifier for the output assets. Must not contain whitespace. Recommended: use "-" instead of "_" if needed.')
    key_params.add_argument('--title', type=str, help='(Parameters for --run-cartload2) The title of the output assets. If it contains spaces, quote it, e.g., "Example Title"')
    key_params.add_argument('--desc', type=str, help='(Parameters for --run-cartload2) The description of output assets. If it contains spaces, quote it, e.g., "Example short description"')
    # cell and boundary assets
    key_params.add_argument('--cell-id', type=str, default="xeniumranger", help='(Parameters for --import-cells) Identifier of the cell factor. This will be used as the ID and filename prefix for cell and boundary assets. Must not contain whitespace. (default: xeniumranger).')
    key_params.add_argument('--cell-name', type=str, help='(Parameters for --import-cells) Name of the cell factor. This will be used in the cell assets (defaults to --cell-id)')
    key_params.add_argument('--tsv-cmap', type=str, default=None, help='(Parameters for --import-cells) Path to the TSV file containing the color map for the cell clusters')
    # images
    key_params.add_argument('--image-ids', type=str, default=["DAPI_OME", "BOUNDARY_OME", "INTERIOR_RNA_OME", "INTERIOR_PROTEIN_OME", "DAPI_MIP_OME"], nargs="+", help='(Parameters for --import-images) One or more image IDs to be used in the output PMTiles. (default: "DAPI_OME", "BOUNDARY_OME", "INTERIOR_RNA_OME", "INTERIOR_PROTEIN_OME", "DAPI_MIP_OME").')
    key_params.add_argument('--image-colors', type=str, nargs='*', default=[], help='(Parameters for --import-images) One or more colors to be used for OME TIFFs. The order of the colors should match the order of the image IDs in --image-ids. If not set, a default list of colors will be used (palette covers up to 10 images). For more than 10 images, please provide --image-colors explicitly.')
    key_params.add_argument('--all-images', action='store_true', help='If set, deploy all detected images based on the --xenium-ranger-assets file. If set, ignore the --image-ids.')
    key_params.add_argument("--transparent-below", type=int, default=1, help='Set pixels below this value to transparent. The threshold should be between 0 and 255 (default: 1).')
    # aws 
    key_params.add_argument("--s3-bucket", help="(Parameters for --upload-aws; required) AWS S3 bucket (e.g., cartostore).")
    # zenodo
    key_params.add_argument('--zenodo-token', type=str,  help='(Parameters for --upload-zenodo; required) Path to a file containing your Zenodo access token.')
    key_params.add_argument('--zenodo-deposition-id', type=str, default=None,  help='(Parameters for --upload-zenodo; optional) ID of an existing Zenodo deposition to upload files into. If not provided, a new deposition will be created. If the deposition is published, a new version will be automatically created.')
    key_params.add_argument('--zenodo-title', type=str, default=None, help='(Parameters for --upload-zenodo) Title of the deposition. If it contains spaces, quote it, e.g., "Example Zenodo Title". Recommended: provide when creating a new deposition (defaults to --title)') 
    key_params.add_argument('--creators', type=str, nargs='+', default=[], help='(Parameters for --upload-zenodo) List of creators, each enclosed in double quotes, in "Lastname, Firstname" format. Recommended: provide when creating a new deposition')  

    aux_inout_params = parser.add_argument_group("Auxiliary Input/Output File/Directories Parameters", "Manually specify input file names and output directories. If the input file locations (--csv-* and --ome-tifs) are not provided, load from --xenium-ranger-assets. Similarly, if output directories (--sge-dir, --fic-dir, --cart-dir) are not specified, defaults under --out-dir will be used.")
    aux_inout_params.add_argument('--xenium-ranger-assets', type=str, default=None, help='Path to a JSON file containing Xenium Ranger asset paths. This file is written when --load-xenium-ranger is enabled, and is read when no manual input files (--csv-*, --ome-tifs) are specified. Default: <out_dir>/xenium_ranger_assets.json if --out-dir is set.')
    aux_inout_params.add_argument('--csv-transcript', type=str, default=None, help='(Parameters for --sge-convert) Location of the CSV or parquet file in the --xenium-ranger-dir. This transcript file should contain at least coordinates, feature names, and expression counts of the transcripts')
    aux_inout_params.add_argument('--csv-cells', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV or Parquet file containing cell locations in the --xenium-ranger-dir') 
    aux_inout_params.add_argument('--csv-boundaries', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV file containing cell boundary coordinates in the --xenium-ranger-dir')
    aux_inout_params.add_argument('--csv-clust', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV file containing cell cluster assignments in the --xenium-ranger-dir') 
    aux_inout_params.add_argument('--csv-diffexp', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV file with differential expression results in the --xenium-ranger-dir') 
    aux_inout_params.add_argument('--ome-tifs', type=str, nargs="+", default=[], help="(Parameters for --import-images) List of locations of one or more input OME-TIFF(s) in the --xenium-ranger-dir. The order must match the order of the image IDs in --image-ids")
    aux_inout_params.add_argument('--catalog-yaml', type=str, nargs="+", default=[], help="(Parameters for --import-images) List of locations of one or more input OME-TIFF(s) in the --xenium-ranger-dir. The order must match the order of the image IDs in --image-ids")

    # Default settings: * tools will use the PATH: sort, gzip, python3
    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--pmtiles', type=str, help='Path to pmtiles binary from go-pmtiles (default: pmtiles)')
    env_params.add_argument('--gdal_translate', type=str, help='Path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdaladdo', type=str, help='Path to gdaladdo binary (default: gdaladdo)')
    env_params.add_argument('--gdalwarp', type=str, help='Path to gdalwarp binary (default: gdalwarp)')
    env_params.add_argument('--tippecanoe', type=str, help='Path to tippecanoe binary (default: <cartloader_dir>/submodules/tippecanoe/tippecanoe)')
    env_params.add_argument('--spatula', type=str, help='Path to spatula binary (default: spatula)')
    env_params.add_argument('--parquet-tools', type=str, dest='parquet_tools', help='Path to parquet-tools binary (default: parquet-tools)')
    env_params.add_argument('--ficture2', type=str, default=os.path.join(repo_dir, "submodules", "punkst"),  help='Path to punkst(ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--aws', type=str, default="aws", help='The path to aws (default: aws)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    # Required-when validations with helpful CLI errors
    # * Directory/layout requirement
    if not args.out_dir and not (args.sge_dir and args.fic_dir and args.cart_dir):
        parser.error("Either --out-dir must be specified, or all of --sge-dir, --fic-dir, and --cart-dir must be specified.")

    if (not args.xenium_ranger_assets) and args.out_dir:
        args.xenium_ranger_assets = os.path.join(args.out_dir, "xenium_ranger_assets.json")

    # * Manual inputs vs asset JSON are mutually exclusive when --load-xenium-ranger
    manual_mode = any([ args.csv_transcript, args.csv_cells, args.csv_boundaries, args.csv_clust, args.csv_diffexp ] + (args.ome_tifs or []))
    if manual_mode:
        if args.load_xenium_ranger:
            parser.error("Cannot enable both auto-detection model (--load-xenium-ranger) and manual input model (--csv-*, --ome-tifs). Choose one mode")
        if args.all_images:
            parser.error("--all-images cannot be used in manual input model (specifying --ome-tifs)")
        if not args.xenium_ranger_dir:
            parser.error("--xenium-ranger-dir is required in manual input model (specifying --csv-*, --ome-tifs)")
        if not os.path.exists(args.xenium_ranger_dir):
            parser.error(f"Directory not found: {args.xenium_ranger_dir} (--xenium-ranger-dir)")
    else:
        assert args.xenium_ranger_assets, "Provide --xenium-ranger-assets (or --out-dir) when using auto-detection mode"

    # * load xenium_ranger/manual input
    if args.load_xenium_ranger:
        if not args.xenium_ranger_dir:
            parser.error("--xenium-ranger-dir is required when --load-xenium-ranger is set")
        if not os.path.exists(args.xenium_ranger_dir):
            parser.error(f"Directory not found: {args.xenium_ranger_dir} (--xenium-ranger-dir)")
    
    # * general args
    validate_general_args(parser, args)

    return args

def run_xenium(_args):
    # parse argument
    args=parse_arguments(_args)

    # Directories
    os.makedirs(args.out_dir, exist_ok=True)

    sge_dir = os.path.join(args.out_dir, "sge")
    fic_dir = os.path.join(args.out_dir, "ficture2")
    cart_dir = os.path.join(args.out_dir, "cartload2")

    if args.import_cells or args.import_images or args.run_cartload2 or args.import_ext_ficture2:
        os.makedirs(cart_dir, exist_ok=True)
    
    # Manual vs JSON inputs
    manual_inputs = [args.csv_transcript, args.csv_cells, args.csv_boundaries, args.csv_clust, args.csv_diffexp] + (args.ome_tifs or [])
    use_json = not any(x is not None for x in manual_inputs)

    # Common paths
    ranger_dir = args.xenium_ranger_dir
    ranger_assets = args.xenium_ranger_assets
    sge_assets = os.path.join(sge_dir, "sge_assets.json")
    sge_flag = os.path.join(sge_dir, "sge_density_filtering.done" if args.filter_by_density else "sge_convert.done")
    cell_assets = os.path.join(cart_dir, f"{args.cell_id}_assets.json")
    background_assets = []
    background_pmtiles = []
    catalog_yaml = os.path.join(cart_dir, "catalog.yaml")

    # Actions in order
    # * Steps that will be executed directly: load_xenium_ranger, import_xenium_cell, upload_zenodo
    # * Steps that will be executed via make file: sge_convert, run_ficture2, import-images, run_cartload2, upload_aws, 
    if args.load_xenium_ranger:
        print("="*10, flush=True)
        print("Executing --load-xenium-ranger (execute directly)", flush=True)
        print("="*10, flush=True)

        print(f" * Input Ranger Directory: {ranger_dir}", flush=True)
        print(f" * Output JSON: {ranger_assets}", flush=True)
        
        load_xenium_cmd=" ".join([
                "cartloader", "load_xenium_ranger",
                f"--in-dir {ranger_dir}",
                f"--out-json {ranger_assets}",
                f"--unzip-dir {args.out_dir}/unzip" 
            ])

        if os.path.exists(ranger_assets) and not args.restart:
            print(f" * Skip --load-xenium-ranger since the Xenium Ranger raw input assets file ({ranger_assets}) already exists. You can use --restart to force execution of this step.", flush=True)
            print(load_xenium_cmd, flush=True)
        else:
            run_command_w_preq(load_xenium_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    
    if args.sge_convert:
        print("="*10, flush=True)
        print("Executing --sge-convert (execute via make)", flush=True)
        print("="*10, flush=True)

        os.makedirs(sge_dir, exist_ok=True) 
        
        if use_json:
            prereq = [ranger_assets]
        else:
            assert args.csv_transcript, "--csv-transcript is required when using manual inputs for --sge-convert"
            csv_transcript_path=os.path.join(ranger_dir, args.csv_transcript)
            assert os.path.exists(csv_transcript_path), f"File not found: {csv_transcript_path} (provided by --xenium-ranger-dir and --csv-transcript)"
            prereq = [csv_transcript_path]
        
        print(f" * Input JSON: {ranger_assets}" if use_json else f"Input Transcript: {os.path.join(ranger_dir, args.csv_transcript)}", flush=True)
        print(f" * SGE directory: {sge_dir}", flush=True)

        sge_convert_cmd = " ".join([
            "cartloader", "sge_convert",
            f"--platform 10x_xenium",
            f"--out-dir {sge_dir}",
            f"--in-json {ranger_assets}" if use_json else "",
            f"--in-csv {os.path.join(ranger_dir, args.csv_transcript)}" if not use_json else "",
            f"--units-per-um {args.units_per_um}" if args.units_per_um else "",
            f"--filter-by-density" if args.filter_by_density else "",
            f"--exclude-feature-regex {args.exclude_feature_regex}" if args.exclude_feature_regex else "",
            "--sge-visual --north-up", # always north up
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
            ])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["spatula", "parquet_tools", "gdalwarp"])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["restart","n_jobs","threads"])
        
        run_command_w_preq(sge_convert_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)
        
    if args.run_ficture2:
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets]
        stage_run_ficture2(fic_dir, sge_assets, args, prereq)

    if args.import_cells:
        print("="*10, flush=True)
        print("Executing --import-cells (execute directly)", flush=True)
        print("="*10, flush=True)
        
        if use_json:
            prereq = [ranger_assets]
        else:
            prereq = [os.path.join(ranger_dir, x) for x in [args.csv_cells, args.csv_boundaries, args.csv_clust, args.csv_diffexp] if x is not None]

        print(f" * Input JSON: {ranger_assets}" if use_json else f"Input Directory: {ranger_dir}")    

        import_cell_cmd=" ".join([
            "cartloader", "import_xenium_cell",
            # actions
            "--all",
            f"--update-catalog" if not args.run_cartload2 and os.path.exists(catalog_yaml) else "",
            # output
            f"--outprefix {cart_dir}/{args.cell_id}",
            # (conditional) input
            f"--in-json {ranger_assets}" if use_json else "",
            f"--in-dir {ranger_dir}" if not use_json else "",
            f"--csv-cells {args.csv_cells}" if (not use_json) and args.csv_cells else "",
            f"--csv-boundaries {args.csv_boundaries}" if (not use_json) and args.csv_boundaries else "",
            f"--csv-clust {args.csv_clust}" if (not use_json) and args.csv_clust else "",
            f"--csv-diffexp {args.csv_diffexp}" if (not use_json) and args.csv_diffexp else "",
            # key params
            f"--id {args.cell_id}",
            f"--name {args.cell_name}" if args.cell_name else "",
            f"--tsv-cmap {args.tsv_cmap}" if args.tsv_cmap else "",
            f"--tippecanoe {args.tippecanoe}" if args.tippecanoe else "",
            f"--threads {args.threads}" if args.threads else "",
        ])

        if os.path.exists(cell_assets) and not args.restart:
            print(f" * Skip --import-cells since the Xenium Ranger cell assets file ({cell_assets}) already exists. You can use --restart to force execution of this step.", flush=True)
            print(import_cell_cmd, flush=True)
        else:
            run_command_w_preq(import_cell_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    
    if args.import_images:  ## TBC: currently, only support OME tifs
        print("="*10, flush=True)
        print("Executing --import-images (execute via make)", flush=True)
        print("="*10, flush=True)

        args.image_ids = validate_imageid_args(args, ranger_assets)
        args.image_colors = validate_imagecol_args(args)
        
        # (xenium-specific) ome_tifs when not use_json
        tif_paths=[]
        if not use_json:
            assert args.ome_tifs, "--ome-tifs is required when --import-images is set with manual input mode. Otherwise, set --load-xenium-ranger to automatically detect the images."
            assert_unique(args.ome_tifs, "--ome-tifs")
            # length
            if len(args.ome_tifs) < len(args.image_ids):
                raise ValueError(f"--ome-tifs ({len(args.ome_tifs)}) must be >= --image-ids ({len(args.image_ids)}); ensure paths align in order")
            # all items in OME-TIFFs should exist:
            for tif_loc in args.ome_tifs:
                tif_path = os.path.join(ranger_dir, tif_loc)
                assert os.path.exists(tif_path), f"File not found: {tif_path} (provided by --xenium-ranger-dir and --ome-tifs)"
                tif_paths.append(tif_path)
        
        # cmd and execution
        stage_import_images(cart_dir, args, ranger_assets, use_json, tif_paths, catalog_yaml)

        background_assets  = [ os.path.join(cart_dir, f"{image_id}_assets.json") for image_id in args.image_ids]
        background_pmtiles = [ os.path.join(cart_dir, f"{image_id}.pmtiles")     for image_id in args.image_ids]

    if args.run_cartload2:   
        # prerequisites
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets] 
        
        # cmd and execution
        stage_run_cartload2(cart_dir, fic_dir, cell_assets, background_assets, args, prereq)
        
    if args.upload_aws:
        # prerequisites
        prereq = [catalog_yaml]
        # cmd and execution
        stage_upload_aws(cart_dir, args, prereq=prereq)

    if args.upload_zenodo:
        # prerequisites
        prereq = [catalog_yaml]
        # Define flag files for skip logic
        zenodo_cartload_flag = f"{cart_dir}/cartload.zenodo.done"
        zenodo_basemap_flags = [f"{bg}.zenodo.done" for bg in background_pmtiles]  if background_pmtiles else []
        flag = zenodo_basemap_flags + [zenodo_cartload_flag]
        # cmd and execution
        stage_upload_zenodo(cart_dir, args, prereq=prereq, flag=flag)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
