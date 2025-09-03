import sys, os, argparse, logging, subprocess, inspect
import pandas as pd
from pathlib import Path
import hashlib

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import add_param_to_cmd, cmd_separator, run_command_w_preq, load_file_to_dict, assert_unique
from cartloader.utils.pipeline_helper import validate_general_args, validate_imageid_args, validate_imagecol_args, stage_run_ficture2, stage_run_cartload2, stage_upload_aws, stage_upload_zenodo, stage_import_images

def parse_arguments(_args):
    """
    Process Visium HD output
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Process 10x Visium HD output from Space Ranger, and return rastered and tiled sources for CartoScope")

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run mode. Print all commands without executing them.')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    cmd_params = parser.add_argument_group("Commands")
    cmd_params.add_argument('--load-space-ranger', action='store_true', default=False, help='Automatically detect Space Ranger output files in --space-ranger-dir and save their paths to --space-ranger-assets (JSON).')
    cmd_params.add_argument('--sge-convert', action='store_true', default=False, help='Convert SGE files from VisiumHD MEX format into the cartloader-standardized format.')
    cmd_params.add_argument('--run-ficture2', action='store_true', default=False, help='Run FICTURE analysis on standardized SGE using the latest punkst version.')
    cmd_params.add_argument('--import-ext-ficture2', action='store_true', default=False, help='Import a set of existing FICTURE results.')
    cmd_params.add_argument('--import-cells', action='store_true', default=False, help='Import cell-based analysis results from Space Ranger, and add the asset information to catalog.yaml if exists.')
    cmd_params.add_argument('--import-images', action='store_true', default=False, help='Import images as background (e.g., DAPI) and add asset information to catalog.yaml if present.')
    cmd_params.add_argument('--run-cartload2', action='store_true', default=False, help='Package SGE into raster tiles (and import FICTURE results if --run-ficture2 is enabled). A catalog.yaml file will be created to summarize all assets.')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload output assets into AWS')
    cmd_params.add_argument('--upload-zenodo', action='store_true', default=False, help='Upload output assets into Zenodo')

    inout_params = parser.add_argument_group("Input/Output Parameters", 'cartloader determines input files in two ways:	1. Manual input mode (explicit) — If any argument in "Auxiliary Input/Output File/Directories" are specified with --xenium-ranger-dir; 2. Auto-detection mode (implicit) — If no manual input is given, cartloader reads input file paths from --xenium-ranger-assets, which is created by --load-xenium-ranger.')
    inout_params.add_argument('--space-ranger-dir', type=str, help='(Required if manual input mode is enabled or running --load-space-ranger)  Path to the Space Ranger output directory containing transcript, cell, boundary, cluster, and image files.')
    inout_params.add_argument('--out-dir', type=str, help='Output root directory. The output for --sge-convert is written to <out_dir>/sge; for --run-ficture2 to <out_dir>/ficture2; and for --run-cartload2 to <out_dir>/cartload2.')
    inout_params.add_argument('--ext-fic-dir', type=str, help='(Required if --import-ext-ficture2 is enabled) Path to an external FICTURE directory for loading external FICTURE assets.')

    key_params = parser.add_argument_group("Key Parameters")
    key_params.add_argument('--space-ranger-assets', type=str, default=None, help='Path to a JSON file containing VisiumHD Space Ranger paths. This file is written when --load-space-ranger is enabled, and is read when no manual input files (--csv-*, --btf-tifs) are specified (default: <out_dir>/space_ranger_assets.json).')
    # sge convert
    key_params.add_argument('--units-per-um', type=float, default=None, help='(Parameters for --sge-convert) Coordinate unit per um in the input files (default: 1.00).')  
    key_params.add_argument('--filter-by-density', action='store_true', default=False, help='(Parameters for --sge-convert) Filter SGE during format conversion by density (default: False).')
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
    key_params.add_argument('--cell-id', type=str, default="spaceranger", help='(Parameters for --import-cells) Identifier of the cell factor. This will be used as the ID and filename prefix for cell and boundary assets. Must not contain whitespace. (default: spaceranger).')
    key_params.add_argument('--cell-name', type=str, help='(Parameters for --import-cells) Name of the cell factor. This will be used in the cell assets (defaults to --cell-id)')
    key_params.add_argument('--tsv-cmap', type=str, default=None, help='(Parameters for --import-cells) Path to the TSV file containing the color map for the cell clusters')
    # images
    key_params.add_argument('--image-ids', type=str, default=["DAPI_OME", "BOUNDARY_OME", "INTERIOR_RNA_OME", "INTERIOR_PROTEIN_OME", "DAPI_MIP_OME"], nargs="+", help='(Parameters for --import-images) One or more image IDs to be used in the output PMTiles. (default: "DAPI_OME", "BOUNDARY_OME", "INTERIOR_RNA_OME", "INTERIOR_PROTEIN_OME", "DAPI_MIP_OME").')
    key_params.add_argument('--image-colors', type=str, nargs='*', default=[], help='(Parameters for --import-images) One or more colors to be used for OME TIFFs. The order of the colors should match the order of the image IDs in --image-ids. If not set, a default list of colors will be used (palette covers up to 10 images). For more than 10 images, please provide --image-colors explicitly.')
    key_params.add_argument('--all-images', action='store_true', help='If set, deploy all detected images based on the --space-ranger-assets file. If set, ignore the --image-ids.')
    key_params.add_argument("--transparent-below", type=int, default=1, help='Set pixels below this value to transparent. The threshold should be between 0 and 255 (default: 1).')
    # aws 
    key_params.add_argument("--s3-bucket", help="(Parameters for --upload-aws; required) AWS S3 bucket (e.g., cartostore).")
    # zenodo
    key_params.add_argument('--zenodo-token', type=str,  help='(Parameters for --upload-zenodo; required) Path to a file containing your Zenodo access token.')
    key_params.add_argument('--zenodo-deposition-id', type=str, default=None,  help='(Parameters for --upload-zenodo; optional) ID of an existing Zenodo deposition to upload files into. If not provided, a new deposition will be created. If the deposition is published, a new version will be automatically created.')
    key_params.add_argument('--zenodo-title', type=str, default=None, help='(Parameters for --upload-zenodo) Title of the deposition. If it contains spaces, quote it, e.g., "Example Zenodo Title". Recommended: provide when creating a new deposition (defaults to --title)') 
    key_params.add_argument('--creators', type=str, nargs='+', default=[], help='(Parameters for --upload-zenodo) List of creators, each enclosed in double quotes, in "Lastname, Firstname" format. Recommended: provide when creating a new deposition')  

    aux_inout_params = parser.add_argument_group("Manual Input Files Parameters", "Manual input file locations (paths relative to --space-ranger-dir). If any of the parameters are defined, manual input mode is enabled.")
    aux_inout_params.add_argument('--mex-transcript', type=str, default=None, help='(Parameters for --sge-convert) Location of the transcript file in CSV or parquet format in the --space-ranger-dir. This transcript file should contain at least coordinates, feature names, and expression counts of the transcripts')
    aux_inout_params.add_argument('--json-scale', type=str, default=None, help='(Parameters for --sge-convert) Location of the JSON file containing spatial scale factor(s) in the --space-ranger-dir.')
    aux_inout_params.add_argument('--parquet-position', type=str, default=None, help='(Parameters for --sge-convert) Location of the Parquet file containing the barcode coordinates in the --space-ranger-dir.')
    aux_inout_params.add_argument('--geojson-cells', type=str, default=None, help='(Parameters for --import-cells) Location of the GEOJSON file containing cell locations in the --space-ranger-dir')
    aux_inout_params.add_argument('--mtx-cells', type=str, default=None, help='(Parameters for --import-cells) Location of the cell feature matrix in MTX format')
    aux_inout_params.add_argument('--csv-clust', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV file containing cell cluster assignments in the --space-ranger-dir') 
    aux_inout_params.add_argument('--csv-diffexp', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV file with differential expression results in the --space-ranger-dir') 
    aux_inout_params.add_argument('--btf-tifs', type=str, nargs="+", default=[], help="(Parameters for --import-images) One or more BTF TIFF input paths in the --space-ranger-dir. The order should match the image IDs in --image-ids.")

    # Default settings: tools will use the PATH: sort, gzip, python3
    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--pmtiles', type=str, help='Path to pmtiles binary from go-pmtiles (default: pmtiles)')
    env_params.add_argument('--gdal_translate', type=str, help='Path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdaladdo', type=str, help='Path to gdaladdo binary (default: gdaladdo)')
    env_params.add_argument('--gdalwarp', type=str, help='Path to gdalwarp binary (default: gdalwarp)')
    env_params.add_argument('--tippecanoe', type=str, help='Path to tippecanoe binary (default: <cartloader_dir>/submodules/tippecanoe/tippecanoe)')
    env_params.add_argument('--spatula', type=str, help='Path to spatula binary (default: spatula)')
    env_params.add_argument('--parquet-tools', type=str, dest='parquet_tools', help='Path to parquet-tools binary (default: parquet-tools)')
    env_params.add_argument('--ficture2', type=str, default=os.path.join(repo_dir, "submodules", "punkst"),  help='Path to punkst (ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--aws', type=str, default="aws", help='The path to aws (default: aws)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args(_args)

    # Required-when validations with helpful CLI errors
    # * Directory/layout requirement
    if not args.out_dir:
        parser.error("--out-dir must be specified.")

    if not args.space_ranger_assets:
        args.space_ranger_assets = os.path.join(args.out_dir, "space_ranger_assets.json")

    # * Manual inputs vs asset JSON are mutually exclusive when --load-space-ranger
    manual_mode = any([ args.mex_transcript, args.json_scale, args.geojson_cells, args.mtx_cells, args.csv_clust, args.csv_diffexp] + (args.btf_tifs or []))
    if manual_mode:
        if args.load_space_ranger:
            parser.error("Cannot enable both auto-detection mode (--load-space-ranger) and manual input mode (--mex-transcript, --json-scale, --geojson-cells, --mtx-cells, --csv-*, --btf-tifs). Choose one mode.")
        if args.all_images:
            parser.error("--all-images cannot be used in manual input mode (when specifying --btf-tifs).")
        if not args.space_ranger_dir:
            parser.error("--space-ranger-dir is required in manual input mode (when specifying --mex-transcript, --json-scale, --geojson-cells, --mtx-cells, --csv-*, --btf-tifs).")
        if not os.path.exists(args.space_ranger_dir):
            parser.error(f"Directory not found: {args.space_ranger_dir} (--space-ranger-dir)")
    else:
        assert args.space_ranger_assets, "Provide --space-ranger-assets (or --out-dir) when using auto-detection mode"

    # * load space_ranger/manual input
    if args.load_space_ranger:
        if not args.space_ranger_dir:
            parser.error("--space-ranger-dir is required when --load-space-ranger is set")
        if not os.path.exists(args.space_ranger_dir):
            parser.error(f"Directory not found: {args.space_ranger_dir} (--space-ranger-dir)")
    
    # * general args
    validate_general_args(parser, args)

    return args

def run_visiumhd(_args):
    # parse argument
    args=parse_arguments(_args)
    
    # Directories
    os.makedirs(args.out_dir, exist_ok=True)

    sge_dir = os.path.join(args.out_dir, "sge")
    fic_dir = os.path.join(args.out_dir, "ficture2")
    cart_dir = os.path.join(args.out_dir, "cartload2")
    
    if args.import_cells or args.import_images or args.run_cartload2 or args.import_ext_ficture2:
        os.makedirs(cart_dir, exist_ok=True)
    
    # use json or not
    manual_inputs = [ args.mex_transcript, args.json_scale, args.geojson_cells, args.mtx_cells, args.csv_clust, args.csv_diffexp] + (args.btf_tifs or [])
    use_json = not any(x is not None for x in manual_inputs)

    # Common paths
    ranger_dir = args.space_ranger_dir
    ranger_assets = args.space_ranger_assets
    sge_assets = os.path.join(sge_dir, "sge_assets.json")
    sge_flag = os.path.join(sge_dir, "sge_density_filtering.done" if args.filter_by_density else "sge_convert.done")
    cell_assets = os.path.join(cart_dir, f"{args.cell_id}_assets.json")
    background_assets = []
    background_pmtiles = []
    catalog_yaml = os.path.join(cart_dir, "catalog.yaml")

    # Actions in order
    # * Steps that will be executed directly: load_space_ranger, import_*_cell, upload_zenodo
    # * Steps that will be executed via make file: sge_convert, run_ficture2, import-images, run_cartload2, upload_aws, 
    if args.load_space_ranger:
        print("="*10, flush=True)
        print("Executing --load-space-ranger (execute directly)", flush=True)
        print("="*10, flush=True)

        print(f" * Input Ranger directory: {ranger_dir}", flush=True)
        print(f" * Output JSON: {ranger_assets}", flush=True)
        
        load_space_cmd=" ".join([
                "cartloader", "load_space_ranger",
                f"--in-dir {args.space_ranger_dir}",
                f"--out-json {args.space_ranger_assets}",
                f"--unzip-dir {args.out_dir}/unzip" 
            ])
    
        if os.path.exists(ranger_assets) and not args.restart:
            print(f" * Skip --load-space-ranger since the Space Ranger raw input assets file ({args.space_ranger_assets}) already exists. You can use --restart to force execution of this step.", flush=True)
            print(load_space_cmd, flush=True)
        else:
            run_command_w_preq(load_space_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
        
    sge_assets = os.path.join(sge_dir, "sge_assets.json")
    if args.sge_convert:
        print("="*10, flush=True)
        print("Executing --sge-convert (execute via make)", flush=True)
        print("="*10, flush=True)

        sge_flag = f"{sge_dir}/sge_convert.done" if not args.filter_by_density else f"{sge_dir}/sge_density_filtering.done"

        if use_json:
            prereq = [args.space_ranger_assets]
        else:
            sge_convert_input={
                "--mex-transcript": args.mex_transcript,
                "--parquet-position": args.parquet_position
            }

            if args.json_scale:
                sge_convert_input.update({"--json-scale": args.json_scale})
            
            if not args.json_scale and not args.units_per_um:
                raise ValueError(
                        "Cannot find scale factor. Specify it using one of the following options:\n"
                        "  1) Enable automatic detection with --load-space-ranger\n"
                        "  2) Specify the file path with --json-scale\n"
                        "  3) Provide the scale factor directly with --units-per-um"
                    )

            prereq=[]
            for key, value in sge_convert_input:
                assert value, f"{key} is required when using manual inputs for --sge-convert"
                file_path=os.path.join(ranger_dir, value)
                assert os.path.exists(file_path), f"File not found: {file_path} ({key})"
                prereq.append(file_path)
            
        print(f" * SGE directory: {sge_dir}", flush=True)
        print(f" * Input JSON: {args.space_ranger_assets}" if use_json else f"Input transcript directory: {os.path.join(ranger_dir, args.mex_transcript)}\nInput position parquet: {os.path.join(ranger_dir, args.parquet_position)}\nInput scale: {os.path.join(ranger_dir, args.json_scale) if args.json_scale else 'N/A'}", flush=True)
        
        sge_convert_cmd = " ".join([
            "cartloader", "sge_convert",
            f"--platform 10x_visium_hd",
            f"--out-dir {sge_dir}",
            f"--in-json {args.space_ranger_assets}" if use_json else "",
            f"--in-mex {os.path.join(ranger_dir, args.mex_transcript)}" if not use_json else "",
            f"--pos-parquet {os.path.join(ranger_dir, args.parquet_position)}" if not use_json else "",
            f"--scale-json {os.path.join(ranger_dir, args.json_scale)}" if (not use_json and args.json_scale) else "",
            f"--units-per-um {args.units_per_um}" if (not use_json and args.units_per_um) else "",
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
            prereq = [args.space_ranger_assets]
        else:
            prereq = [os.path.join(args.space_ranger_dir, x) for x in [args.geojson_cells, args.mtx_cells ,args.csv_clust, args.csv_diffexp] if x is not None]

        print(f" * Input JSON: {args.space_ranger_assets}" if use_json else f"Input directory: {args.space_ranger_dir}")        
        import_cell_cmd=" ".join([
            "cartloader", "import_visiumhd_cell",
            # actions
            "--all",
            f"--update-catalog" if not args.run_cartload2 and os.path.exists(catalog_yaml) else "",
            # output
            f"--outprefix {cart_dir}/{args.cell_id}",
            # (conditional) input
            f"--in-json {args.space_ranger_assets}" if use_json else "",
            f"--in-dir {args.space_ranger_dir}" if not use_json else "",
            f"--geojson-cells {args.geojson_cells}" if (not use_json) and args.geojson_cells else "",
            f"--mtx-cells {args.mtx_cells}" if (not use_json) and args.mtx_cells else "",
            f"--csv-clust {args.csv_clust}" if (not use_json) and args.csv_clust else "",
            f"--csv-diffexp {args.csv_diffexp}" if (not use_json) and args.csv_diffexp else "",
            # key params
            f"--id {args.cell_id}" if args.cell_id else "",
            f"--name {args.cell_name}" if args.cell_name else "",
            f"--tsv-cmap {args.tsv_cmap}" if args.tsv_cmap else "",
            f"--tippecanoe {args.tippecanoe}" if args.tippecanoe else "",
            f"--threads {args.threads}" if args.threads else "",
        ])

        if os.path.exists(cell_assets) and not args.restart:
            print(f" * Skip --import-cells since the Space Ranger cell assets file ({cell_assets}) already exists. You can use --restart to force execution of this step.", flush=True)
            print(import_cell_cmd, flush=True)
        else:
            run_command_w_preq(import_cell_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    
    if args.import_images:  ## TBC: currently supports only OME-TIFFs
        print("="*10, flush=True)
        print("Executing --import-images (execute via make)", flush=True)
        print("="*10, flush=True)

        args.image_ids = validate_imageid_args(args, ranger_assets)
        args.image_colors = validate_imagecol_args(args)
        
        # (visiumhd-specific) btf when not use_json
        tif_paths=[]
        if not use_json:
            assert args.btf_tifs, "--btf-tifs is required when --import-images is set with manual input mode. Otherwise, set --load-space-ranger to automatically detect the images."
            assert_unique(args.btf_tifs, "--btf-tifs")
            # length
            if len(args.btf_tifs) < len(args.image_ids):
                raise ValueError(f"--btf-tifs ({len(args.btf_tifs)}) must be >= --image-ids ({len(args.image_ids)}); ensure paths align in order")
            # all items in BTF TIFFs should exist:
            for tif_loc in args.btf_tifs:
                tif_path = os.path.join(ranger_dir, tif_loc)
                assert os.path.exists(tif_path), f"File not found: {tif_path} (provided by --space-ranger-dir and --btf-tifs)"
                tif_paths.append(tif_path)

        stage_import_images(cart_dir, args, ranger_assets, use_json, tif_paths, catalog_yaml)

        background_assets  = [ os.path.join(cart_dir, f"{image_id}_assets.json") for image_id in args.image_ids]
        background_pmtiles = [ os.path.join(cart_dir, f"{image_id}.pmtiles")     for image_id in args.image_ids]

    if args.run_cartload2:   
        # prerequisites
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets] 
        stage_run_cartload2(cart_dir, fic_dir, cell_assets, background_assets, args, prereq)
        
    if args.upload_aws:
        # prerequisites
        prereq = [catalog_yaml]

        stage_upload_aws(cart_dir, args, prereq=prereq)

    if args.upload_zenodo:
        # prerequisites
        prereq = [catalog_yaml]

        # Define flag files for skip logic
        zenodo_cartload_flag = f"{cart_dir}/cartload.zenodo.done"
        zenodo_basemap_flags = [f"{bg}.zenodo.done" for bg in background_pmtiles]  if background_pmtiles else []
        flag = zenodo_basemap_flags + [zenodo_cartload_flag]

        stage_upload_zenodo(cart_dir, args, prereq=prereq, flag=flag)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
