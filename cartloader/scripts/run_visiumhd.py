import sys, os, argparse, logging, subprocess, inspect
import pandas as pd
from pathlib import Path
import hashlib

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import add_param_to_cmd, cmd_separator, run_command_w_preq, load_file_to_dict
from cartloader.scripts.run_xenium import stage_run_ficture2, stage_run_cartload2, stage_upload_aws, stage_upload_zenodo

def parse_arguments(_args):
    """
    Process Xenium output
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Process 10X Xenium output from Xenium Ranger, and return rastered and tiled sources for CartoScope")

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run mode. Print all commands without executing them.')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    cmd_params = parser.add_argument_group("Commands")
    cmd_params.add_argument('--load-space-ranger', action='store_true', default=False, help='Automatically detect Xenium Ranger output files in --space-ranger-dir and save their paths to --space-ranger-assets (JSON).')
    cmd_params.add_argument('--sge-convert', action='store_true', default=False, help='Convert SGE files from Xenium format into the cartloader-standardized format.')
    cmd_params.add_argument('--run-ficture2', action='store_true', default=False, help='Run FICTURE analysis on standardized SGE using the latest punkst version.')
    cmd_params.add_argument('--import-ext-ficture2', action='store_true', default=False, help='Import an existing FICTURE results.')
    cmd_params.add_argument('--import-cells', action='store_true', default=False, help='Import the Xenium Ranger output, including cells, boundaries, clusters, and differentially expression genes. If catalog.yaml exists in --cart-dir, update it with asset information.')
    cmd_params.add_argument('--import-images', action='store_true', default=False, help='Import images as background, such as DAPI, and add the asset information to catalog.yaml in --cart-dir if the yaml file exists.')
    cmd_params.add_argument('--run-cartload2', action='store_true', default=False, help='Package SGE into raster tiles (and import FICTURE results if --run-ficture2 is enabled). A catalog.yaml file will be created to summarize all assets.')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload output assets into AWS')
    cmd_params.add_argument('--upload-zenodo', action='store_true', default=False, help='Upload output assets into Zenodo')

    inout_params = parser.add_argument_group("Input/Output Parameters", 'Cartloader determines input files in two ways:	1. Manual input (explicit) — If any --csv-* or --ome-tifs arguments in "Auxiliary Input/Output File/Directories" are specified with --space-ranger-dir, these file paths relative to are used directly. 2. Asset JSON (implicit) — If no manual input is given, cartloader reads input file paths from --space-ranger-assets, which is created by --load-space-ranger with --space-ranger-dir.')
    inout_params.add_argument('--space-ranger-dir', type=str, help='(Required if using --load-space-ranger is enabled or any manual input by --csv-* or --ome-tifs) Path to the Xenium Ranger output directory containing transcript, cell, boundary, cluster, and image files.')
    inout_params.add_argument('--out-dir', type=str, help='Path to output root directory. If set, --sge-dir, --fic-dir, --cart-dir, --space-ranger-assets default to <out_dir>/sge, <out_dir>/ficture2, <out_dir>/cartload2, and <out_dir>/space_ranger_assets.json. To use a custom directory layout, skip --out-dir and define those paths manually.')
    inout_params.add_argument('--ext-fic-dir', type=str, help='(Required if --import-ext-ficture2 is enabled) Path to an external FICTURE directory for loading external FICTURE assets.')

    key_params = parser.add_argument_group("Key Parameters")
    key_params.add_argument('--space-ranger-assets', type=str, default=None, help='Path to a JSON file containing VisiumHD Space Ranger paths. This file is written when --load-space-ranger is enabled, and is read when no manual input files (--csv-*, --ome-tifs) are specified. Default: <out_dir>/space_ranger_assets.json if --out-dir is set.')
    # sge convert
    key_params.add_argument('--units-per-um', type=float, default=None, help='(Parameters for --sge-convert) Coordinate unit per um in the input files (default: 1.00).')  
    key_params.add_argument('--filter-by-density', action='store_true', default=False, help='(Parameters for --sge-convert) Filter SGE from format conversion by density (default: False).')
    key_params.add_argument('--exclude-feature-regex', type=str, default=None, help='(Parameters for --sge-convert) A regex pattern of feature/gene names to be excluded (default: "^(BLANK_|DeprecatedCodeword_|NegCon|UnassignedCodeword_)")')
    # ficture2 parameters
    key_params.add_argument('--width', type=str, default=None, help='(Parameters for --run-ficture2; required) Comma-separated hexagon flat-to-flat widths (in um) for LDA training and projection default: None)') ## Same width will be used for both train width and projection width
    key_params.add_argument('--n-factor', type=str, default=None, help='(Parameters for --run-ficture2; required) Comma-separated list of factor counts for LDA training.')
    key_params.add_argument('--colname-feature', type=str, default='gene', help='(Parameters for --run-ficture2 or --run-cartload2) Input/output Column name for gene name (default: gene)')
    key_params.add_argument('--colname-count', type=str, default='count', help='(Parameters for --run-ficture2 or --run-cartload2) Column name for feature counts (default: count)')
    # run_cartload2 parameters
    key_params.add_argument('--id', type=str, help='(Parameters for --run-cartload2; required) Identifier for the output assets. Must not contain whitespace. Recommended: use "-" instead of "_" if needed.')
    key_params.add_argument('--title', type=str, help='(Parameters for --run-cartload2) The title of the output assets')
    key_params.add_argument('--desc', type=str, help='(Parameters for --run-cartload2) The description of output assets')
    # cell and boundary assets
    key_params.add_argument('--cell-id', type=str, default="xeniumranger", help='(Parameters for --import-cells) Identifier of the cell factor. This will be used as the ID and filename prefix for cell and boundary assets (default: space_ranger).')
    key_params.add_argument('--cell-name', type=str, help='(Parameters for --import-cells) Name of the cell factor. This will be used in the cell assets (defaults to --cell-id)')
    key_params.add_argument('--tsv-cmap', type=str, default=None, help='(Parameters for --import-cells) Path to the TSV file containing the color map for the cell clusters')
    # images
    key_params.add_argument('--image-ids', type=str, default=["BTF"], nargs="+", help='(Parameters for --import-images) One or more image IDs to be used in the output PMTiles. (default: "BTF").')
    key_params.add_argument('--image-colors', type=str, default=[], help='(Parameters for --import-images) One or more colors to be used for OME TIFFs. The order of the colors should match the order of the image IDs in --image-ids. If not set, a default list of colors will be used.')
    key_params.add_argument('--skip-missing-images', action='store_true', help='If set, skip image IDs that do not exist, showing a warning instead of raising an error.')
    key_params.add_argument("--transparent-below", type=int, default=1, help='Set pixels below this value to transparent. The threshold should be between 0 and 255 (default: 1).')

    # aws 
    key_params.add_argument("--s3-bucket", help="(Parameters for --upload-aws; required) AWS S3 bucket (e.g., cartostore).")
    # zenodo
    key_params.add_argument('--zenodo-token', type=str,  help='(Parameters for --upload-zenodo; required) Path to a file containing your Zenodo access token.')
    key_params.add_argument('--zenodo-deposition-id', type=str, default=None,  help='(Parameters for --upload-zenodo; optional) ID of an existing Zenodo deposition to upload files into. If not provided, a new deposition will be created. If the deposition is published, a new version will be automatically created.')
    key_params.add_argument('--zenodo-title', type=str, default=None, help='(Parameters for --upload-zenodo) Title of the deposition. Recommended: provide when creating a new deposition (defaults to --title)') #Required if creating a new deposition or the existing deposition does not have a title
    key_params.add_argument('--creators', type=str, nargs='+', default=[], help='(Parameters for --upload-zenodo) List of creators, each enclosed in double quotes, in "Lastname, Firstname" format. Recommended: provide when creating a new deposition')  #Required if creating a new deposition or the existing deposition do not have creator information.

    aux_inout_params = parser.add_argument_group("Auxiliary Input/Output File/Directories Parameters", "Manually specify input file names and output directories. If the input file locations (--csv-* and --ome-tifs) are not provided, load from --space-ranger-assets. Similarly, if output directories (--sge-dir, --fic-dir, --cart-dir) are not specified, defaults under --out-dir will be used.")
    aux_inout_params.add_argument('--mex-transcript', type=str, default=None, help='(Parameters for --sge-convert) Location of the CSV or parquet file in the --space-ranger-dir. This transcript file should contain at least coordinates, feature names, and expression counts of the transcripts')
    aux_inout_params.add_argument('--json-scale', type=str, default=None, help='(Parameters for --sge-convert) Location of the JSON file containing scale factor information in the --space-ranger-dir') 
    aux_inout_params.add_argument('--geojson-cells', type=str, default=None, help='(Parameters for --import-cells) Location of the GEOJSON file containing cell locations in the --indir')
    aux_inout_params.add_argument('--mtx-cell-ftr', type=str, default=None, help='(Parameters for --import-cells) Location of the cell feature matrix in MTX format')
    aux_inout_params.add_argument('--csv-clust', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV file containing cell cluster assignments in the --space-ranger-dir') 
    aux_inout_params.add_argument('--csv-diffexp', type=str, default=None, help='(Parameters for --import-cells) Location of the CSV file with differential expression results in the --space-ranger-dir') 
    aux_inout_params.add_argument('--ome-tifs', type=str, nargs="+", default=[], help="(Parameters for --import-images) List of locations of one or more input ome tiff(s) in the --space-ranger-dir. The order of the locations should match the order of the image IDs in --image-id")
    aux_inout_params.add_argument('--sge-dir', type=str, help='Path to input/output SGE directory to store or to provide SGE and a YAML/JSON file summarizing the path of transcript, feature, and minmax. Recommended: use --out-dir for consistent directory structure.')
    aux_inout_params.add_argument('--fic-dir', type=str, help='Path to FICTURE directory to store or to provide FICTURE output with a YAML/JSON file file summarizing the SGE and the FICTURE parameters. Recommended: use --out-dir for consistent directory structure.')
    aux_inout_params.add_argument('--cart-dir', type=str, help='Path to cartload directory to store or host the deployed assets, including cell and boundary assets, SGE assets, FICTURE assets, and histology assets. Recommended: use --out-dir for consistent directory structure.')
    
    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    # Default settings:
    # * tools will use the PATH: sort, gzip, python3
    env_params.add_argument('--pmtiles', type=str, help='Path to pmtiles binary from go-pmtiles (default: pmtiles)')
    env_params.add_argument('--gdal_translate', type=str, help='Path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdaladdo', type=str, help='Path to gdaladdo binary (default: gdaladdo)')
    env_params.add_argument('--tippecanoe', type=str, help='Path to tippecanoe binary (default: <cartloader_dir>/submodules/tippecanoe/tippecanoe)')
    env_params.add_argument('--spatula', type=str, help='Path to spatula binary (default: spatula)')
    env_params.add_argument('--ficture2', type=str, default=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "submodules", "punkst"),  help='Path to punkst(ficture2) repository (default: <cartloader_dir>/submodules/punkst)')
    env_params.add_argument('--aws', type=str, default="aws", help='The path to aws (default: aws)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)

def run_visiumhd(_args):
    # parse argument
    args=parse_arguments(_args)
    
    # actions:
    assert not (args.run_ficture2 and args.import_ext_ficture2), "Cannot run both --run-ficture2 and --import-ext-ficture2 at the same time. Please choose one of them."

    # dirs
    # make sure that the user either specified --out-dir or all of --sge-dir, --fic-dir, and --cart-dir
    if not args.out_dir and (not args.sge_dir or not args.fic_dir or not args.cart_dir):
        raise ValueError("Either --out-dir must be specified, or all of --sge-dir, --fic-dir, and --cart-dir must be specified.")
    
    if args.out_dir:
        sge_dir = os.path.join(args.out_dir, "sge") if args.sge_dir is None else args.sge_dir
        fic_dir = os.path.join(args.out_dir, "ficture2") if args.fic_dir is None else args.fic_dir
        cart_dir = os.path.join(args.out_dir, "cartload2") if args.cart_dir is None else args.cart_dir
    else:
        sge_dir = args.sge_dir
        fic_dir = args.fic_dir
        cart_dir = args.cart_dir

    # use json or not
    manual_inputs = [args.mex_transcript, args.geojson_cells, args.mex_cell_ftr, args.csv_clust, args.csv_diffexp] + args.ome_tifs

    if args.load_space_ranger and any(x is not None for x in manual_inputs):
        raise ValueError("Cannot use --load-space-ranger with manually specified inputs (--csv-transcript, --geojson-cells, --csv-boundaries, --csv-clust, or --csv-diffexp). These options are mutually exclusive.")

    use_json = not any(x is not None for x in manual_inputs)

    if (not args.space_ranger_assets) and args.out_dir:
        args.space_ranger_assets = os.path.join(args.out_dir, "space_ranger_assets.json")
    
    # catalog yaml
    catalog_yml = os.path.join(cart_dir, "catalog.yaml")

    # Actions by orders
    # * Steps that will executed directly: load_space_ranger, import_*_cell, upload_zenodo
    # * Steps that will executed via make file: sge_convert, run_ficture2, import-images, run_cartload2, upload_aws, 

    if args.load_space_ranger:
        print("="*10, flush=True)
        print("Executing --load-space-ranger (execute directly)", flush=True)
        print("="*10, flush=True)

        assert args.space_ranger_assets, "Missing path to the VisiumHD Space Ranger Asset JSON file. Specify it using either --out-dir or --space-ranger-assets"
        assert args.space_ranger_dir and os.path.exists(args.space_ranger_dir), "Please specify the directory hosting the Xenium Ranger output files using --space-ranger-dir"

        print(f" * VisiumHD Space Ranger Directory: {args.space_ranger_dir}", flush=True)
        
        load_xenium_cmd=" ".join([
                "cartloader", "detect_visiumhd_output",
                f"--in-dir {args.space_ranger_dir}",
                f"--out-json {args.space_ranger_assets}",
                f"--unzip-dir {args.out_dir}/unzip",
            ])
        if os.path.exists(args.space_ranger_assets) and not args.restart:
            print(f" * Skip --load-space-ranger since the Xenium Ranger raw input assets file {args.space_ranger_assets} already exists. You can use --restart to force execute this step.", flush=True)
            print(load_xenium_cmd, flush=True)
        else:
            run_command_w_preq(load_xenium_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    
    sge_assets = os.path.join(sge_dir, "sge_assets.json")
    if args.sge_convert:
        print("="*10, flush=True)
        print("Executing --sge-convert (execute via make)", flush=True)
        print("="*10, flush=True)

        sge_flag = f"{sge_dir}/sge_convert.done" if not args.filter_by_density else f"{sge_dir}/sge_density_filtering.done"

        if use_json:
            assert args.space_ranger_assets, "Missing path to the VisiumHD Space Ranger Asset JSON file. Specify it using either --out-dir or --space-ranger-assets"
            prereq = [args.space_ranger_assets]
        else:
            assert args.space_ranger_dir and os.path.exists(args.space_ranger_dir), "Please specify the directory hosting the VisiumHD Space Ranger output files using --space-ranger-dir"
            prereq = [os.path.join(args.space_ranger_dir, args.mex_transcript)]
        
        print(f" * SGE directory: {sge_dir}", flush=True)
        print(f" * Input JSON: {args.space_ranger_assets}" if use_json else f"Input MEX: {args.mex_transcript}\nInput Scale:{args.json_scale}", flush=True)
        
        if not use_json and not args.json_scale and not args.units_per_um:
            raise ValueError(
                    "Cannot find scale factor. Specify it using one of the following options:\n"
                    "  1) Enable automatic detection with --load-space-ranger\n"
                    "  2) Specify the file path with --json-scale\n"
                    "  3) Provide the scale factor directly with --units-per-um"
                )
        
        sge_convert_cmd = " ".join([
            "cartloader", "sge_convert",
            f"--platform 10x_vsisium_hd",
            f"--out-dir {sge_dir}",
            f"--in-json {args.space_ranger_assets}" if use_json else "",
            f"--in-mex {args.mex_transcript}" if not use_json else "",
            f"--scale_json {args.json_scale}" if (not use_json and args.json_scale) else "",
            f"--units-per-um {args.units_per_um}" if (not use_json and args.units_per_um) else "",
            f"--filter-by-density" if args.filter_by_density else "",
            f"--exclude-feature-regex {args.exclude_feature_regex}" if args.exclude_feature_regex else "",
            "--sge-visual",
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
            ])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["spatula", "parquet_tools", "gdalwarp"])
        sge_convert_cmd = add_param_to_cmd(sge_convert_cmd, args, ["restart","n_jobs","threads"])
        
        run_command_w_preq(sge_convert_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)
        
    cell_assets = os.path.join(cart_dir, f"{args.cell_id}_assets.json")
    if args.import_cells:
        print("="*10, flush=True)
        print("Executing --import-cells (execute directly)", flush=True)
        print("="*10, flush=True)
        
        # Use cell_assets as the target and flag. import_xenium_cell is executed directly instead of via Makefile
        assert args.cell_id, "Please specify the --cell-id parameter for importing cells and boundaries"
        
        if use_json:
            assert args.space_ranger_assets, "Missing path to the VisiumHD Space Ranger JSON file. Specify it using either --out-dir or --space-ranger-assets"
            prereq = [args.space_ranger_assets]
        else:
            assert args.space_ranger_dir and os.path.exists(args.space_ranger_dir), "Please specify the directory hosting the Xenium Ranger output files using --space-ranger-dir"
            prereq = [os.path.join(args.space_ranger_dir, x) for x in [args.geojson_cells, args.mtx_cell_ftr, args.json_scale ,args.csv_clust, args.csv_diffexp] if x is not None]

        print(f" * Input JSON: {args.space_ranger_assets}" if use_json else f"Input Directory: {args.space_ranger_dir}")        
        import_cell_cmd=" ".join([
            "cartloader", "import_visiumhd_cell",
            # actions
            "--all",
            f"--update-catalog" if not args.run_cartload2 and os.path.exists(catalog_yml) else "",
            # output
            f"--outprefix {cart_dir}/{args.cell_id}",
            # (conditional) input
            f"--in-json {args.space_ranger_assets}" if use_json else "",
            f"--in-dir {args.space_ranger_dir}" if not use_json else "",
            f"--geojson-cells {args.geojson_cells}" if (not use_json) and args.geojson_cells else "",
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
            print(f" * Skip --import-xenium-output since the Xenium Ranger cell assets file ({cell_assets}) already exists. You can use --restart to force execute this step.", flush=True)
            print(import_cell_cmd, flush=True)
        else:
            run_command_w_preq(import_cell_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    
    background_assets = []
    background_pmtiles = []
    if args.import_images:  ## TBC: currently, only support OME tifs
        print("="*10, flush=True)
        print("Executing --import-images (execute via make)", flush=True)
        print("="*10, flush=True)
        
        colors = [     ## A list of colors that work both for dark and light backgrounds
            "#1f77b4",  # Blue
            "#ff7f0e",  # Orange
            "#2ca02c",  # Green
            "#d62728",  # Red
            "#9467bd",  # Purple
            "#8c564b",  # Brown
            "#e377c2",  # Pink
            "#7f7f7f",  # Gray
            "#bcbd22",  # Olive
            "#17becf",  # Cyan
        ]

        if not use_json:
            assert args.ome_tifs, "Please specify the --ome-tifs parameter to import images or set --load-space-ranger to automatically detect the images"
            assert not args.skip_missing_images, "Cannot use (TBC)"

        assert args.image_ids, "Please specify at least one image ID via --image-ids parameter to import images"

        if args.skip_missing_images: # skip_missing_images can only be used with use_json
            if args.dry_run:
                print("* Skip removing unavailable image IDs for --dry-run")
            else:
                # load the assets file
                assert args.space_ranger_assets, "Missing path to the VisiumHD Space Ranger JSON file. Specify it using either --out-dir or --space-ranger-assets"
                print("* Filtering image IDs to only remain available ones (--skip-missing-images)")
                space_ranger_data=load_file_to_dict(args.space_ranger_assets)
                # drop the pairs with None values
                space_ranger_data={k: v for k, v in space_ranger_data.items() if v is not None}
                avail_keys=list(space_ranger_data.keys())
                raw_image_ids=args.image_ids
                args.image_ids = [img_id for img_id in raw_image_ids if img_id in avail_keys]

                print(f"    - Raw image IDs (N={len(raw_image_ids)}): {raw_image_ids}")
                print(f"    - Available image IDs (N={len(args.image_ids)}): {args.image_ids}")
        
        if not args.image_colors:
            args.image_colors = colors[:len(args.image_ids)]
        assert len(args.image_ids) <= len(args.image_colors), "Please specify image colors more or equal to the number of color image IDs"
        
        for i in range(len(args.image_ids)):
            image_id = args.image_ids[i]
            image_color = args.image_colors[i].replace("#","")

            if use_json:
                assert args.space_ranger_assets, "Missing path to the VisiumHD Space Ranger JSON file. Specify it using either --out-dir or --space-ranger-assets"
                prereq = [args.space_ranger_assets]
                ome_path = None  # Not needed, but safe to define for later use
            else:
                assert args.space_ranger_dir and os.path.exists(args.space_ranger_dir), "Please specify the directory hosting the Xenium Ranger output files using --space-ranger-dir"
                ome_path = os.path.join(args.space_ranger_dir, args.ome_tifs[i])
                prereq = [ome_path]
                assert os.path.exists(ome_path), f"Input OME TIFF file {ome_path} does not exist. Please check the path."

            print(f" * Image ID: {image_id}", flush=True)
            print(f" * Image color: {image_color}", flush=True)

            import_image_cmd = " ".join([
                "cartloader", "import_image",
                # actions
                "--ome2png",
                "--png2pmtiles",
                f"--update-catalog" if not args.run_cartload2 and os.path.exists(catalog_yml) else "",
                # input
                f"--in-json {args.space_ranger_assets}" if use_json else "",
                f"--in-img {ome_path}" if not use_json else "",
                # output
                f"--out-dir {cart_dir}",
                # id and color
                f"--img-id {image_id}",
                f"--colorize \"{image_color}\"",
                # # aux params
                # f"--upper-thres-quantile 0.95",
                # f"--lower-thres-quantile 0.5",
                f"--transparent-below {args.transparent_below}" if args.transparent_below else "",
                f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",

            ])
            import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["pmtiles", "gdaladdo"])
            import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["restart","n_jobs"])
            
            run_command_w_preq(import_image_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)

            background_assets.append(os.path.join(cart_dir, f"{image_id}_assets.json"))
            background_pmtiles.append(os.path.join(cart_dir, f"{image_id}.pmtiles"))

    if args.run_ficture2:
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets]
        stage_run_ficture2(fic_dir, sge_assets, args, prereq)
        
    if args.run_cartload2:   
        # prerequisites
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets] 
        stage_run_cartload2(cart_dir, fic_dir, cell_assets, background_assets, args, prereq)
        
    if args.upload_aws:
        prereq = [catalog_yml]
        stage_upload_aws(cart_dir, args, prereq=prereq)

    if args.upload_zenodo:
        prereq = [catalog_yml]

        # Define flag files for skip logic
        zenodo_cartload_flag = f"{cart_dir}/cartload.zenodo.done"
        zenodo_basemap_flags = [f"{bg}.zenodo.done" for bg in background_pmtiles]
        flag = zenodo_basemap_flags + [zenodo_cartload_flag]

        stage_upload_zenodo(cart_dir, args, prereq=prereq, flag=flag)
        

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
