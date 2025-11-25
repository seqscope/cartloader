import sys, os, argparse, logging, subprocess, inspect
import pandas as pd
from pathlib import Path
import hashlib

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import add_param_to_cmd, cmd_separator, run_command_w_preq, load_file_to_dict, assert_unique
from cartloader.utils.pipeline_helper import resolve_image_plan, validate_general_args, validate_imageid_args, validate_imagecol_args, validate_imageloc_args, stage_run_ficture2, stage_run_cartload2, stage_upload_aws, stage_upload_zenodo, stage_import_images, resolve_square_plan, stage_import_squares

manual_str="--mex-transcript, --json-scale, --parquet-position, --geojson-cells, --mtx-cells, --csv-*, --tifs"

def parse_arguments(_args):
    """
    Process Visium HD output
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Process 10x Visium HD output from Space Ranger, and return rastered and tiled sources for CartoScope (load → sge_convert → ficture2 → images/cells → cartload2 → uploads)"
    )

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Print commands without executing')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of parallel jobs to run (default: 1)')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    cmd_params = parser.add_argument_group("Commands")
    cmd_params.add_argument('--load-space-ranger', action='store_true', default=False, help='Detect Space Ranger outputs in --space-ranger-dir and write --space-ranger-assets (JSON)')
    cmd_params.add_argument('--sge-convert', action='store_true', default=False, help='Convert SGE from Visium HD format to cartloader-compatible format')
    cmd_params.add_argument('--run-ficture2', action='store_true', default=False, help='Run FICTURE (punkst) analysis on the converted SGE. Note, --decode-scale is set to 2 by default for Visium HD')
    cmd_params.add_argument('--import-ext-ficture2', action='store_true', default=False, help='Import a set of existing FICTURE results')
    cmd_params.add_argument('--import-cells', action='store_true', default=False, help='Import Space Ranger cell analysis into tiles and update catalog.yaml (if present)')
    cmd_params.add_argument('--import-squares', action='store_true', default=False, help='Import Space Ranger square bins results into tiles and update catalog.yaml (if present)')
    cmd_params.add_argument('--import-images', action='store_true', default=False, help='Import background images (e.g., DAPI) into raster tiles and update catalog.yaml (if present)')
    cmd_params.add_argument('--run-cartload2', action='store_true', default=False, help='Package SGE into raster tiles (import FICTURE results if --run-ficture2). Write catalog.yaml summarizing assets')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload output assets into AWS')
    cmd_params.add_argument('--upload-zenodo', action='store_true', default=False, help='Upload output assets into Zenodo')

    inout_params = parser.add_argument_group(
        "Input/Output Parameters",
        'Two input modes: 1) Manual — set --space-ranger-dir and any "Manual Input Files Parameters" args. 2) Auto-detect — leave them unset; read from --space-ranger-assets (created by --load-space-ranger)'
    )
    inout_params.add_argument('--space-ranger-dir', type=str, help='Path to the Space Ranger output directory containing transcript, cell, boundary, cluster, and image files  (required if manual input mode is enabled or if --load-xenium-ranger)')
    inout_params.add_argument('--out-dir', type=str, required=True, help='Path to output directory. Stores converted SGE, FICTURE results, raster tiles in <out_dir>/sge, <out_dir>/ficture2, and <out_dir>/cartload2')
    inout_params.add_argument('--space-ranger-assets', type=str, default=None, help=f'Path to a JSON file containing Space Ranger asset paths. Written by --load-space-ranger; read when auto-detection mode is on (default: <out_dir>/space_ranger_assets.json)')
    inout_params.add_argument('--ext-fic-dir', type=str, help='Path to an external FICTURE directory for loading external FICTURE assets (required if --import-ext-ficture2)')

    # Parameters for --sge-convert
    sge_params = parser.add_argument_group("Parameters for --sge-convert")
    sge_params.add_argument('--units-per-um', type=float, default=None, help='Coordinate unit per um in raw SGE (default: 1.00)')  
    sge_params.add_argument('--filter-by-density', action='store_true', default=False, help='Enable to filter SGE by density')
    sge_params.add_argument('--exclude-feature-regex', type=str, default=None, help='Regex for feature names to exclude when running --sge-convert(default: "^(BLANK_|DeprecatedCodeword_|NegCon|UnassignedCodeword_)")')

    # Parameters for --run-ficture2
    fic_params = parser.add_argument_group("Parameters for --run-ficture2")
    fic_params.add_argument('--width', type=str, default=None, help='Comma-separated hexagon flat-to-flat widths (in um) for LDA training and projection (required if --run-ficture2)') ## Same width will be used for both train width and projection width
    fic_params.add_argument('--n-factor', type=str, default=None, help='Comma-separated list of factor counts for LDA training (required if --run-ficture2)')
    fic_params.add_argument('--colname-feature', type=str, default='gene', help='Column name for feature name (used with --run-ficture2 and --run-cartload2; default: gene)')
    fic_params.add_argument('--colname-count', type=str, default='count', help='Column name for UMI counts (used with --run-ficture2 and --run-cartload2; default: count)')
    fic_params.add_argument('--fic-include-feature-regex', type=str, default=None, help='Regex of feature names to include when running --run-ficture2')
    fic_params.add_argument('--fic-exclude-feature-regex', type=str, default=None, help='Regex of feature names to exclude when running --run-ficture2, e.g., apply "^(mt-.*$|Gm\\d+$)" for mouse datasets to exclude mitochondrial gene and Pseudogenes')

    # Parameters for --run-cartload2
    cart_params = parser.add_argument_group("Parameters for --run-cartload2")
    cart_params.add_argument('--id', type=str, help='Identifier for output assets; no whitespace; prefer "-" over "_" (required if --run-cartload2)')
    cart_params.add_argument('--title', type=str, help='Asset human-readable title. Quote if contains spaces, e.g., "Example Title"')
    cart_params.add_argument('--desc', type=str, help='Asset description. Quote if contains spaces, e.g., "Example short description"')

    # Parameters for --import-cells
    cells_params = parser.add_argument_group("Parameters for --import-cells")
    cells_params.add_argument('--cell-id', type=str, default="spaceranger", help='Identifier of Space Ranger cell results; no whitespace (default: spaceranger)')
    cells_params.add_argument('--cell-name', type=str, help='Name of Space Ranger cell results (defaults to --cell-id)')
    cells_params.add_argument('--tsv-cmap', type=str, default=None, help='Path to a color map TSV file for cell clusters')

    # Parameters for --import-squares
    squares_params = parser.add_argument_group("Parameters for --import-squares")
    squares_params.add_argument('--square-id', type=str, default="spaceranger", help='Identifier of Space Ranger cell results; no whitespace (default: spaceranger)')
    squares_params.add_argument('--use-parquet-tools', action='store_true', help='Use parquet-tools instead of polars/pigz for parquet to csv conversion (default: False). parquet-tools may be slower for large files.')

    # Parameters for --import-images
    images_params = parser.add_argument_group("Parameters for --import-images", ('Two input modes: 1) use --image-ids select images; 2) use --all-images to deploy all images from --space-ranger-assets. '
                                                                                 'Note, run_visiumhd is currently designed for H&E images in BTF format with OME-XML metadata, i.e., it always define --georef-detect OME. '
                                                                                 'For other images, consider run "cartloader import_image" directly with appropriate parameters.'))
    images_params.add_argument('--image-ids', type=str, default=["HnE"], nargs="+", help='One or more image IDs to import (default: ["HnE"])')
    images_params.add_argument('--all-images', action='store_true', help='Enable to deploy all images from --space-ranger-assets regardless --image-ids')

    # Parameters for --upload-aws or --upload-zenodo
    upload_params = parser.add_argument_group("Parameters for --upload-aws and --upload-zenodo")
    upload_params.add_argument("--s3-dir", help="AWS S3 directory to host output, e.g., s3://cartostore/test (required if --upload-aws)")
    upload_params.add_argument('--zenodo-token', type=str,  help='Path to a file containing your Zenodo access token (required if --upload-zenodo)')
    upload_params.add_argument('--zenodo-deposition-id', type=str, default=None,  help='Existing deposition ID. If published, creates a new version; omit to create a new deposition.')
    upload_params.add_argument('--zenodo-title', type=str, default=None, help='Deposition title; always quote the value, e.g., "Example Zenodo Title" (defaults to --title)') 
    upload_params.add_argument('--creators', type=str, nargs='+', default=[], help='List of creators as "Lastname, Firstname" (each quoted)')  

    aux_inout_params = parser.add_argument_group("Manual Input Files Parameters", "Manually specify input file locations under --space-ranger-dir; providing any of the following arguments enables manual mode.")
    aux_inout_params.add_argument('--mex-transcript', type=str, default=None, help='Location of the CSV or parquet transcript file under --space-ranger-dir containing coordinates, feature names, and counts (used with --sge-convert)')
    aux_inout_params.add_argument('--json-scale', type=str, default=None, help='Location of the JSON with spatial scale factor(s) under --space-ranger-dir (used with --sge-convert)')
    aux_inout_params.add_argument('--parquet-position', type=str, default=None, help='Location of the Parquet file containing barcode coordinates under --space-ranger-dir (used with --sge-convert)')
    aux_inout_params.add_argument('--geojson-cells', type=str, default=None, help='Location of the GEOJSON file containing cell locations under --space-ranger-dir (used with --import-cells)')
    aux_inout_params.add_argument('--mtx-cells', type=str, default=None, help='Location of the cell feature matrix in MTX format under --space-ranger-dir (used with --import-cells)')
    aux_inout_params.add_argument('--csv-clust', type=str, default=None, help='Location of the CSV file containing cell cluster assignments under --space-ranger-dir (used with --import-cells)') 
    aux_inout_params.add_argument('--csv-diffexp', type=str, default=None, help='Location of the CSV file with differential expression results under --space-ranger-dir (used with --import-cells)') 
    aux_inout_params.add_argument('--tifs', type=str, nargs="+", default=[], help="List of locations of one or more input BTF TIFF(s) under --space-ranger-dir. The order must match the order of the image IDs in --image-ids (used with --import-images)")
    aux_inout_params.add_argument('--square-input', type=str, nargs="+", default=[], help="List of locations of one or more input square bins under --space-ranger-dir. For each bin, provide its size and path to each file in the following format: bin_size,bin_pos_parquet,bin_csv_cluster,bin_csv_diffexp,bin_scale_json,bin_csv_umap (used with --import-squares)")

    # Default settings: tools will use the PATH: sort, gzip, python3
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
    if not args.out_dir:
        parser.error("--out-dir must be specified.")

    if not args.space_ranger_assets:
        args.space_ranger_assets = os.path.join(args.out_dir, "space_ranger_assets.json")

    # * Manual inputs vs asset JSON are mutually exclusive when --load-space-ranger
    manual_mode = any([ args.mex_transcript, args.json_scale, args.geojson_cells, args.mtx_cells, args.csv_clust, args.csv_diffexp] + (args.tifs or []))
    
    if manual_mode:
        if args.load_space_ranger:
            parser.error(f"Cannot enable both auto-detection mode (--load-space-ranger) and manual input mode (when specifying any of {manual_str}). Choose one mode.")
        if args.all_images:
            parser.error(f"--all-images cannot be used in manual input mode (when specifying any of {manual_str}).")
        if not args.space_ranger_dir:
            parser.error(f"--space-ranger-dir is required in manual input mode (when specifying any of {manual_str}).")
        if not os.path.exists(args.space_ranger_dir):
            parser.error(f"Directory not found: {args.space_ranger_dir} (--space-ranger-dir)")

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
    manual_inputs = [ args.mex_transcript, args.json_scale, args.geojson_cells, args.mtx_cells, args.csv_clust, args.csv_diffexp] + (args.tifs or [])
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
            print(f"* Skip --load-space-ranger since the Space Ranger raw input assets file ({args.space_ranger_assets}) already exists. You can use --restart to force execution of this step.\n", flush=True)
            # print("\n", flush=True)
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
                        "  1) Enable automatic detection with --load-space-ranger, which will detect scale JSON file\n"
                        "  2) Specify path to the scale JSON file with --json-scale\n"
                        "  3) Provide the scale factor directly with --units-per-um"
                    )

            prereq=[]
            for key, value in sge_convert_input.items():
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
            f"--exclude-feature-regex \"{args.exclude_feature_regex}\"" if args.exclude_feature_regex else "",
            "--sge-visual --north-up", # always north up
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
            f"--use-parquet-tools" if args.use_parquet_tools else "",
            ])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["spatula", "parquet_tools", "gdalwarp"])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["restart","n_jobs"])
        
        run_command_w_preq(sge_convert_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)
        
    if args.run_ficture2:
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets]
        stage_run_ficture2(fic_dir, sge_assets, args, prereq, "--decode-scale 2")

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
            print(f" * Skip --import-cells since the Space Ranger cell assets file ({cell_assets}) already exists. You can use --restart to force execution of this step.\n", flush=True)
            # print("\n", flush=True)
            print(import_cell_cmd, flush=True)
        else:
            run_command_w_preq(import_cell_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    
    if args.import_images:  ## TBC: currently supports only OME-TIFFs
        print("="*10, flush=True)
        print("Executing --import-images (execute via make)", flush=True)
        print("="*10, flush=True)

        image_ids = validate_imageid_args(args.image_ids, ranger_assets, all_images=args.all_images, dry_run=args.dry_run)
        image_colors = []
        
        # (visiumhd-specific) tifs when not use_json
        if not use_json:
            if len(args.tifs) < len(image_ids):
                raise ValueError(f"--tifs ({len(args.tifs)}) must be >= --image-ids ({len(image_ids)}); ensure paths align in order")
            tif_paths = validate_imageloc_args(ranger_dir, tifs_loc=args.tifs, 
                            loc_label="--tifs", ranger_dir_label="--space-ranger-dir")
        else:
            tif_paths=[]

        # image plan (id/color/use_json/img_path/ranger_assets/prereq)
        image_plans = resolve_image_plan(image_ids, image_colors, use_json, ranger_assets, tif_paths)
        
        # cmd and execution
        stage_import_images(cart_dir, args, image_plans, 
                            ome2png=False, 
                            transparent_below=None,
                            georef_detect="OME",
                            update_catalog=True if (not args.run_cartload2 and os.path.exists(catalog_yaml)) else False
                            )

        background_assets  = [ os.path.join(cart_dir, f"{image_id}_assets.json") for image_id in image_ids]
        background_pmtiles = [ os.path.join(cart_dir, f"{image_id}.pmtiles")     for image_id in image_ids]

    if args.import_squares:
        print("="*10, flush=True)
        print("Executing --import-squares (execute directly)", flush=True)
        print("="*10, flush=True)

        square_plans = resolve_square_plan(ranger_dir, use_json, ranger_assets, args.square_input)
        square_plans = sorted(square_plans, key=lambda x: x["bin_size"])
        stage_import_squares(cart_dir, args, square_plans,
                            update_catalog=True if (not args.run_cartload2 and os.path.exists(catalog_yaml)) else False
                            )
        
        square_assets = [ os.path.join(cart_dir, f"{args.square_id}-sq{spec["bin_size"]:03d}_assets.json") for spec in square_plans]

    if args.run_cartload2:   
        # prerequisites
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets] 

        

        stage_run_cartload2(cart_dir, fic_dir, sge_dir, 
                            cell_assets, square_assets, background_assets, 
                            args, 
                            prereq)
        
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
