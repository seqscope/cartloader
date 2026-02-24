import sys, os, argparse, logging, subprocess, inspect, re
from pathlib import Path
import hashlib
from typing import Iterable, List, Optional

# from cartloader.utils.minimake import minimake
from cartloader.utils.utils import add_param_to_cmd, cmd_separator, run_command_w_preq, load_file_to_dict, assert_unique
from cartloader.utils.pipeline_helper import resolve_image_plan, validate_general_args, validate_imageid_args, validate_imagecol_args, validate_imageloc_args, stage_run_ficture2, stage_run_cartload2, stage_upload_aws, stage_upload_zenodo, stage_import_images


def parse_arguments(_args):
    """
    Process Xenium output
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Process 10x Xenium output from Xenium Ranger, and return rastered and tiled sources for CartoScope (load → sge_convert → ficture2 → images/cells → cartload2 → uploads")

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Print commands without executing')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of parallel jobs to run (default: 1)')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')
    
    cmd_params = parser.add_argument_group("Commands")
    cmd_params.add_argument('--load-xenium-ranger', action='store_true', default=False, help='Detect Xenium Ranger outputs in --xenium-ranger-dir and write --xenium-ranger-assets (JSON)')
    cmd_params.add_argument('--sge-convert', action='store_true', default=False, help='Convert SGE from Xenium format to cartloader-compatible format')
    cmd_params.add_argument('--run-ficture2', action='store_true', default=False, help='Run FICTURE (punkst) analysis on the converted SGE')
    cmd_params.add_argument('--import-ext-ficture2', action='store_true', default=False, help='Import a set of existing FICTURE results')
    cmd_params.add_argument('--import-cells', action='store_true', default=False, help='Import Xenium Ranger cell analysis into tiles and update catalog.yaml (if present)')
    cmd_params.add_argument('--import-images', action='store_true', default=False, help='Import background images (e.g., DAPI) into raster tiles and update catalog.yaml (if present)')
    cmd_params.add_argument('--run-cartload2', action='store_true', default=False, help='Package SGE into raster tiles (import FICTURE results if --run-ficture2). Write catalog.yaml summarizing assets')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload output assets into AWS')
    cmd_params.add_argument('--upload-zenodo', action='store_true', default=False, help='Upload output assets into Zenodo')

    inout_params = parser.add_argument_group("Input/Output Parameters", 'Two input modes: 1) Manual — set --xenium-ranger-dir and any "Manual Input Parameters" args. 2) Auto-detect — leave them unset; read from --xenium-ranger-assets (created by --load-xenium-ranger)')
    inout_params.add_argument('--xenium-ranger-dir', type=str, help='Path to the Xenium Ranger output directory containing transcript, cell, boundary, cluster, and image files (required if manual input mode is enabled or if --load-xenium-ranger)')
    inout_params.add_argument('--out-dir', type=str, required=True, help='Path to output directory. Stores converted SGE, FICTURE results, raster tiles in <out_dir>/sge, <out_dir>/ficture2, and <out_dir>/cartload2')
    inout_params.add_argument('--xenium-ranger-assets', type=str, default=None, help='Path to a JSON file containing Xenium Ranger outputs paths. Written by --load-space-ranger; read when auto-detection mode is on (default: <out_dir>/xenium_ranger_assets.json)')
    inout_params.add_argument('--ext-fic-dir', type=str, help='Path to an external FICTURE directory for loading external FICTURE assets (required if --import-ext-ficture2)')

    # Key Parameters (split by command)
    sge_params = parser.add_argument_group("Parameters for --sge-convert")
    sge_params.add_argument('--units-per-um', type=float, default=None, help='Coordinate unit per um in raw SGE (default: 1.00)')  
    sge_params.add_argument('--filter-by-density', action='store_true', default=False, help='Enable to filter SGE by density')
    sge_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex for feature names to exclude (default: "^(BLANK_|DeprecatedCodeword_|NegCon|UnassignedCodeword_)")')

    fic_params = parser.add_argument_group("Parameters for --run-ficture2")
    fic_params.add_argument('--width', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (in um) for LDA training and projection (required if --run-ficture2)') ## Same width will be used for both train width and projection width
    fic_params.add_argument('--n-factor', type=str, default=None, help='Comma-separated list of factor counts for LDA training (required if --run-ficture2)')
    fic_params.add_argument('--colname-feature', type=str, default='gene', help='Column name for feature name (used with --run-ficture2 and --run-cartload2; default: gene)')
    fic_params.add_argument('--colname-count', type=str, default='count', help='Column name for UMI counts (used with --run-ficture2 and --run-cartload2; default: count)')
    fic_params.add_argument('--fic-include-feature-regex', type=str, default=None, help='Regex of feature names to include when running --run-ficture2')
    fic_params.add_argument('--fic-exclude-feature-regex', type=str, default=None, help='Regex of feature names to exclude when running --run-ficture2, e.g., apply "^(mt-.*$|Gm\\d+$)" for mouse datasets to exclude mitochondrial gene and Pseudogenes')

    cart_params = parser.add_argument_group("Parameters for --run-cartload2")
    cart_params.add_argument('--id', type=str, help='Identifier for output assets; no whitespace; prefer "-" over "_" (required if --run-cartload2)')
    cart_params.add_argument('--title', type=str, help='Asset human-readable title. Quote if contains spaces, e.g., "Example Title"')
    cart_params.add_argument('--desc', type=str, help='Asset description. Quote if contains spaces, e.g., "Example short description"')

    cells_params = parser.add_argument_group("Parameters for --import-cells")
    cells_params.add_argument('--cell-id', type=str, default="xeniumranger", help='Identifier of Xenium Ranger cell results; no whitespace (default: xeniumranger)') # This will be used as the ID and filename prefix for cell and boundary assets. 
    cells_params.add_argument('--cell-name', type=str, help='Name of Xenium Ranger cell results (defaults to --cell-id)')
    cells_params.add_argument('--tsv-cmap', type=str, default=None, help='Path to a color map TSV file for cell clusters')

    images_params = parser.add_argument_group("Parameters for --import-images", "Two input modes: 1) use --image-ids select images; 2) use --all-images to deploy all images from --xenium-ranger-assets")
    images_params.add_argument('--image-ids', type=str, default=["DAPI_OME", "BOUNDARY_OME", "INTERIOR_RNA_OME", "INTERIOR_PROTEIN_OME", "DAPI_MIP_OME"], nargs="+", help='One or more image IDs to import (default: "DAPI_OME", "BOUNDARY_OME", "INTERIOR_RNA_OME", "INTERIOR_PROTEIN_OME", "DAPI_MIP_OME").')
    images_params.add_argument('--image-colors', type=str, nargs='*', default=[], help='List of 6-digits HEX RGB codes (e.g., #1f77b4 or 1f77b4). Order matches --image-ids. Defaults cover up to 10 images; provide more if needed.')
    images_params.add_argument('--all-images', action='store_true', help='Enable to deploy all images from --xenium-ranger-assets regardless --image-ids')
    images_params.add_argument("--transparent-below", type=int, default=1, help='Set pixels below this value to transparent (range: 0~255; default: 1)')

    upload_params = parser.add_argument_group("Parameters for --upload-aws and --upload-zenodo")
    upload_params.add_argument("--s3-dir", help="AWS S3 directory to host output, e.g., s3://cartostore/test (required if --upload-aws)")
    upload_params.add_argument("--profile", type=str, default=None, help="AWS CLI profile to use when uploading to S3")
    upload_params.add_argument('--zenodo-token', type=str,  help='Path to a file containing your Zenodo access token (required if --upload-zenodo)')
    upload_params.add_argument('--zenodo-deposition-id', type=str, default=None, help='Existing deposition ID. If published, creates a new version; omit to create a new deposition.')
    upload_params.add_argument('--zenodo-title', type=str, default=None, help='Deposition title; always quote the value, e.g., "Example Zenodo Title" (defaults to --title)') 
    upload_params.add_argument('--creators', type=str, nargs='+', default=[], help='List of creators as "Lastname, Firstname" (each quoted)')  

    aux_inout_params = parser.add_argument_group("Manual Input Parameters", "Manually specify input file locations under --xenium-ranger-dir; providing any of the following arguments enables manual mode.")
    aux_inout_params.add_argument('--csv-transcript', type=str, default=None, help='Location of the CSV or parquet file under --xenium-ranger-dir. This transcript file should contain at least coordinates, feature names, and expression counts (used with --sge-convert)')
    aux_inout_params.add_argument('--csv-cells', type=str, default=None, help='Location of the CSV or Parquet file containing cell locations under --xenium-ranger-dir (used with --import-cells)') 
    aux_inout_params.add_argument('--csv-boundaries', type=str, default=None, help='Location of the CSV file containing cell boundary coordinates under --xenium-ranger-dir (used with --import-cells)')
    aux_inout_params.add_argument('--csv-clust', type=str, default=None, help='Location of the CSV file containing cell cluster assignments under --xenium-ranger-dir (used with --import-cells)') 
    aux_inout_params.add_argument('--csv-diffexp', type=str, default=None, help='Location of the CSV file with differential expression results under --xenium-ranger-dir (used with --import-cells)') 
    aux_inout_params.add_argument('--ome-tifs', type=str, nargs="+", default=[], help="List of locations of one or more input OME-TIFF(s) under --xenium-ranger-dir. The order must match the order of the image IDs in --image-ids (used with --import-images)")

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
    env_params.add_argument('--aws', type=str, default="aws", help='Path to aws CLI (default: aws)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    # Required-when validations with helpful CLI errors
    # * Directory/layout requirement
    if (not args.xenium_ranger_assets) and args.out_dir:
        args.xenium_ranger_assets = os.path.join(args.out_dir, "xenium_ranger_assets.json")

    # * Manual inputs vs auto-detect
    manual_mode = any([ args.csv_transcript, args.csv_cells, args.csv_boundaries, args.csv_clust, args.csv_diffexp ] + (args.ome_tifs or []))
    
    if manual_mode and args.load_xenium_ranger:
        parser.error("Cannot enable both auto-detection model (--load-xenium-ranger) and manual input model (--csv-*, --ome-tifs). Choose one mode")
    
    if manual_mode:            
        if args.all_images:
            parser.error("--all-images cannot be used in manual input model (specifying --ome-tifs)")
        if not args.xenium_ranger_dir:
            parser.error("--xenium-ranger-dir is required in manual input model (specifying --csv-*, --ome-tifs)")
        if not os.path.exists(args.xenium_ranger_dir):
            parser.error(f"Directory not found: {args.xenium_ranger_dir} (--xenium-ranger-dir)")
    
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
            print("\n", flush=True)
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
            f"--exclude-feature-regex \"{args.exclude_feature_regex}\"" if args.exclude_feature_regex else "",
            "--sge-visual --north-up", # always north up
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
            ])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["spatula", "parquet_tools", "gdalwarp"])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["restart","n_jobs"])
        
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
            f"--spatula {args.spatula}" if args.spatula else "",
        ])

        if os.path.exists(cell_assets) and not args.restart:
            print(f" * Skip --import-cells since the Xenium Ranger cell assets file ({cell_assets}) already exists. You can use --restart to force execution of this step.", flush=True)
            print("\n", flush=True)
            print(import_cell_cmd, flush=True)
        else:
            run_command_w_preq(import_cell_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    
    if args.import_images:  ## TBC: currently, only support OME tifs
        print("="*10, flush=True)
        print("Executing --import-images (execute via make)", flush=True)
        print("="*10, flush=True)

        image_ids = validate_imageid_args(args.image_ids, ranger_assets, all_images=args.all_images, dry_run=args.dry_run)
        image_colors = validate_imagecol_args(image_ids, args.image_colors)
        
        # (xenium-specific) ome_tifs when not use_json
        if not use_json:
            if len(args.ome_tifs) < len(image_ids):
                raise ValueError(f"--ome-tifs ({len(args.ome_tifs)}) must be >= --image-ids ({len(image_ids)}); ensure paths align in order")
            tif_paths = validate_imageloc_args(ranger_dir, tifs_loc=args.ome_tifs, 
                            loc_label="--ome-tifs", ranger_dir_label="--xenium-ranger-dir")
        else:
            tif_paths=[]

        # image plan (id/color/use_json/img_path/ranger_assets/prereq)
        image_plans = resolve_image_plan(image_ids, image_colors, use_json, ranger_assets, tif_paths)
        
        # cmd and execution
        stage_import_images(cart_dir, args, image_plans, 
                            ome2png=True, 
                            transparent_below=args.transparent_below,
                            georef_detect=None,
                            update_catalog=True if not args.run_cartload2 and os.path.exists(catalog_yaml) else False
                            )

        background_assets  = [ os.path.join(cart_dir, f"{image_id}_assets.json") for image_id in image_ids]
        background_pmtiles = [ os.path.join(cart_dir, f"{image_id}.pmtiles")     for image_id in image_ids]

    if args.run_cartload2:   
        # prerequisites
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets] 
        
        # cmd and execution
        stage_run_cartload2(cart_dir, fic_dir, sge_dir, cell_assets, [], background_assets, args, prereq)
        
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
