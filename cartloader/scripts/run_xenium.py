import sys, os, argparse, logging, subprocess, inspect
import pandas as pd
from pathlib import Path
import hashlib

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import add_param_to_cmd, cmd_separator, run_command_w_preq


def parse_arguments(_args):
    """
    Process Xenium output
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Process 10X xenium output from Xenium Ranger, and return rastered and tiled sources for CartoScope")

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run mode. If enabled, prints all commands instead of only those that would be executed. Not triggered via the Makefile.')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--threads', type=int, default=4, help='Maximum number of threads per job (for tippecanoe)')

    cmd_params = parser.add_argument_group("Commands")
    cmd_params.add_argument('--load-xenium-ranger', action='store_true', default=False, help='(Shortcut) Automatically detect and summarize the xenium output files from a directory (--xenium-ranger-dir) into a json file. If enabled, no need to manually define the input files using "Auxiliary Input File Name Parameters"')
    cmd_params.add_argument('--sge-convert', action='store_true', default=False, help='Convert SGE files from Xenium format into the cartloader-standardized format.')
    cmd_params.add_argument('--run-ficture2', action='store_true', default=False, help='Run FICTURE analysis on standardized SGE using the latest punkst version.')
    cmd_params.add_argument('--import-ext-ficture2', action='store_true', default=False, help='Import an existing FICTURE results.')
    cmd_params.add_argument('--import-cells', action='store_true', default=False, help='Import the Xenium Ranger output, including cells, boundaries, clusters, and differentially expression genes, and add the asset information to catalog.yaml in --cart-dir if the yaml file exists.')
    cmd_params.add_argument('--import-images', action='store_true', default=False, help='Import images as background, such as DAPI, and add the asset information to catalog.yaml in --cart-dir if the yaml file exists.')
    cmd_params.add_argument('--run-cartload2', action='store_true', default=False, help='Package SGE into raster tiles (and import FICTURE results if --run-ficture2 is enabled). A catalog.yaml file will be created to summarize all assets.')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload output assets into AWS')
    cmd_params.add_argument('--upload-zenodo', action='store_true', default=False, help='Upload output assets into Zenodo')

    inout_params = parser.add_argument_group("Input/Output Parameters", "")
    inout_params.add_argument('--xenium-ranger-dir', type=str, help='(Required if --load-xenium-ranger is enabled or input files are manually provided via --csv-* or --ome-tif) Path to the Xenium Ranger output directory containing transcript, cell, boundary, cluster, and image files.')
    inout_params.add_argument('--root-dir', type=str, help='(Shortcut) Path to a root directory. If set, the paths for --sge-dir, --fic-dir, and --cart-dir default to <root_dir>/sge, <root_dir>/ficture2, and <root_dir>/cartload2.')
    inout_params.add_argument('--xenium-ranger-assets', type=str, default=None, help='Path to a YAML/JSON file to write or read Xenium Ranger asset paths. If --load-xenium-ranger is enabled and this is not provided, defaults to <root_dir>/xenium_ranger_assets.json. Otherwise, specify its path to enable us')
    inout_params.add_argument('--sge-dir', type=str, help='Path to input/output SGE directory to store or to provide SGE and a YAML/JSON file summarizing the path of transcript, feature, and minmax. The file name of this YAML/JSON file is provided via --sge-assets (defaults to <root_dir>/sge).')
    inout_params.add_argument('--fic-dir', type=str, help='Path to FICTURE directory to store or to provide FICTURE output with a YAML/JSON file file summarizing the SGE and the FICTURE parameters. The file name of this YAML/JSON file is specified by --fic-params (defaults to <root_dir>/ficture2).')
    inout_params.add_argument('--cart-dir', type=str, help='Path to cartload directory to store or host the deployed assets, including cell and boundary assets, SGE assets, FICTURE assets, and histology assets. The file name of this YAML/JSON file is specified by --cart-assets (defaults to <root_dir>/cartload2).')
    inout_params.add_argument('--ext-fic-dir', type=str, help='(Required if --) Path to an external FICTURE directory for loading external FICTURE assets.')

    key_params = parser.add_argument_group("Key Parameters")
    # sge convert
    key_params.add_argument('--units-per-um', type=float, default=None, help='(Parameters for --sge-convert) Coordinate unit per um in the input files (default: 1.00).')  
    key_params.add_argument('--filter-by-density', action='store_true', default=False, help='(Parameters for --sge-convert) Filter SGE from format conversion by density (default: False).')
    key_params.add_argument('--exclude-feature-regex', type=str, default=None, help='(Parameters for --sge-convert) A regex pattern of feature/gene names to be excluded. By default, it removes  (default: None)')
    # ficture2 parameters
    key_params.add_argument('--width', type=str, default=None, help='(Parameters for --run-ficture2) Comma-separated hexagon flat-to-flat widths (in um) for LDA training and projection default: None)') ## Same width will be used for both train width and projection width
    key_params.add_argument('--n-factor', type=str, default=None, help='(Parameters for --run-ficture2) Comma-separated list of factor counts for LDA training.')
    key_params.add_argument('--colname-feature', type=str, default='gene', help='(Parameters for --run-ficture2 or --run-cartload2) Input/output Column name for gene name (default: gene)')
    key_params.add_argument('--colname-count', type=str, default='count', help='(Parameters for --run-ficture2 or --run-cartload2) Column name for feature counts')
    # run_cartload2 parameters
    key_params.add_argument('--id', type=str, help='(Parameters for --run-cartload2) The identifier of the output assets. This should has no whitespace')
    key_params.add_argument('--title', type=str, help='(Parameters for --run-cartload2) The title of the output assets')
    key_params.add_argument('--desc', type=str, help='(Parameters for --run-cartload2) The description of output assets')
    # cell and boundary assets
    key_params.add_argument('--cell-id', type=str, default="xeniumranger", help='(Parameters for --import-cells) Identifier of the cell factor. This will be used as the output file name stem for the cell and boundary assets, as well as the ID in the cell assets (default: xenium_ranger).')
    key_params.add_argument('--cell-name', type=str, help='(Parameters for --import-cells) Name of the cell factor. This will be used in the cell assets (defaults to --cell-id)')
    key_params.add_argument('--tsv-cmap', type=str, default=None, help='(Parameters for --import-cells) Path to the TSV file containing the color map for the cell clusters')
    # images
    key_params.add_argument('--image-ids', type=str, default=["DAPI"], nargs="+", help='(Parameters for --import-images) One or more image IDs to be used in the output PMTiles (default: "DAPI").')
    key_params.add_argument('--image-colors', type=str, default=[], help='(Parameters for --import-images) One or more colors to be used for OME TIFFs. The order of the colors should match the order of the image IDs in --image-ids. If not set, a default list of colors will be used. ')
    # aws 
    key_params.add_argument("--s3-bucket", help="(Parameters for --upload-aws) AWS S3 bucket (e.g., cartostore).")
    # zenodo
    key_params.add_argument('--zenodo-token', type=str,  help='(Parameters for --upload-zenodo) Path to a file containing your Zenodo access token.')
    key_params.add_argument('--zenodo-deposition-id', type=str, default=None,  help='(Parameters for --upload-zenodo; optional) ID of an existing Zenodo deposition to upload files into. '
                                                                                        'If the deposition is published, a new version will be automatically created. '
                                                                                        'If not provided, a new deposition will be created.')
    key_params.add_argument('--zenodo-title', type=str, default=None, help='(Parameters for --upload-zenodo) Title of the deposition. Required if creating a new deposition or the existing deposition does not have a title (defaults to --id)')
    key_params.add_argument('--creators', type=str, nargs='+', default=[], help='(Parameters for --upload-zenodo) List of creators in "Lastname, Firstname" format. Required if creating a new deposition or the existing deposition do not have creator information.')

    aux_xenium_params = parser.add_argument_group("Auxiliary Input File Name Parameters", "Used only if neither --load-xenium-ranger is enabled nor an existing assets file is provided via --xenium-ranger-assets.")
    aux_xenium_params.add_argument('--csv-transcript', type=str, default="transcripts.csv.gz", help='(Parameters for --sge-convert) Location of the CSV or parquet file in the --xenium-ranger-dir. This transcript file should contain at least coordinates, feature names, and expression counts of the transcripts')
    aux_xenium_params.add_argument('--csv-cells', type=str, default="cells.csv.gz", help='(Parameters for --import-cells) Location of the CSV or Parquet file containing cell locations in the --xenium-ranger-dir (default: cells.csv.gz)')
    aux_xenium_params.add_argument('--csv-boundaries', type=str, default="cell_boundaries.csv.gz", help='(Parameters for --import-cells) Location of the CSV file containing cell boundary coordinates in the --xenium-ranger-dir (default: cell_boundaries.csv.gz)')
    aux_xenium_params.add_argument('--csv-clust', type=str, default="analysis/clustering/gene_expression_graphclust/clusters.csv", help='(Parameters for --import-cells) Location of the CSV file containing cell cluster assignments in the --xenium-ranger-dir (default: analysis/clustering/gene_expression_graphclust/clusters.csv)')
    aux_xenium_params.add_argument('--csv-diffexp', type=str, default="analysis/diffexp/gene_expression_graphclust/differential_expression.csv", help='(Parameters for --import-cells) Location of the CSV file with differential expression results in the --xenium-ranger-dir (default: analysis/diffexp/gene_expression_graphclust/differential_expression.csv)')
    aux_xenium_params.add_argument('--ome-tifs', type=str, nargs="+", default=[], help="(Parameters for --import-images) List of locations of input ome tiff(s) in the --xenium-ranger-dir.")

    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    # Default settings:
    # * tools will use the PATH: sort, gzip, python3
    env_params.add_argument('--pmtiles', type=str, help='Path to pmtiles binary from go-pmtiles')
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

def run_xenium(_args):
    # parse argument
    args=parse_arguments(_args)
    
    # actions:
    assert not (args.run_ficture2 and args.import_ext_ficture2), "Cannot run both --run-ficture2 and --import-ext-ficture2 at the same time. Please choose one of them."

    # dirs
    # make sure that the user either specified --root-dir or all of --sge-dir, --fic-dir, and --cart-dir
    if not args.root_dir and (not args.sge_dir or not args.fic_dir or not args.cart_dir):
        raise ValueError("Either --root-dir must be specified, or all of --sge-dir, --fic-dir, and --cart-dir must be specified.")
    
    if args.root_dir:
        sge_dir = os.path.join(args.root_dir, "sge") if args.sge_dir is None else args.sge_dir
        fic_dir = os.path.join(args.root_dir, "ficture2") if args.fic_dir is None else args.fic_dir
        cart_dir = os.path.join(args.root_dir, "cartload2") if args.cart_dir is None else args.cart_dir
    else:
        sge_dir = args.sge_dir
        fic_dir = args.fic_dir
        cart_dir = args.cart_dir
    
    # input
    if args.load_xenium_ranger or args.xenium_ranger_assets:
        use_json=True
        if not args.xenium_ranger_assets:
            assert args.root_dir, "The --xenium-ranger-assets must be specified if --load-xenium-ranger is enabled. Please specify the --root-dir or --xenium-ranger-assets."
            args.xenium_ranger_assets = os.path.join(args.root_dir, "xenium_ranger_assets.json")
    else:
        use_json=False
    
    # catalog yaml
    catalog_yml = os.path.join(cart_dir, "catalog.yaml")

    # Actions by orders
    # * Steps that will executed directly: load_xenium_ranger, import_xenium_output, upload_zenodo
    # * Steps that will executed via make file: sge_convert, run_ficture2, import-images, run_cartload2, upload_aws, 

    if args.load_xenium_ranger:
        print("="*10, flush=True)
        print("Executing --load-xenium-ranger (execute directly)", flush=True)
        print("="*10, flush=True)

        assert args.xenium_ranger_assets, "Missing --xenium-ranger-assets. Specify the --xenium-ranger-dir or --root-dir to provide its path."
        print(f" * Xenium Ranger Directory: {args.xenium_ranger_dir}", flush=True)
        load_xenium_cmd=" ".join([
                "cartloader", "detect_xenium_output",
                f"--in-dir {args.xenium_ranger_dir}",
                f"--out-json {args.xenium_ranger_assets}",
                f"--unzip-dir {args.root_dir}/unzip",
            ])
        if os.path.exists(args.xenium_ranger_assets) and not args.restart:
            print(f" * Skip --load-xenium-ranger since the Xenium Ranger raw input assets file {args.xenium_ranger_assets} already exists. You can use --restart to force execute this step.", flush=True)
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
            prereq = [args.xenium_ranger_assets]
        else:
            prereq = [os.path.join(args.xenium_ranger_dir, args.csv_transcript)]
        
        print(f" * SGE directory: {sge_dir}", flush=True)
        print(f" * Input JSON: {args.xenium_ranger_assets}" if use_json else f"Input CSV: {args.csv_transcript}", flush=True)
        sge_convert_cmd = " ".join([
            "cartloader", "sge_convert",
            f"--platform 10x_xenium",
            f"--out-dir {sge_dir}",
            f"--in-json {args.xenium_ranger_assets}" if use_json else "",
            f"--in-csv {args.csv_transcript}" if not use_json else "",
            f"--units-per-um {args.units_per_um}" if args.units_per_um else "",
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
        
        # Use cell_assets as the target and flag. import_xenium_output is executed directly instead of via Makefile
        assert args.cell_id, "Please specify the --cell-id parameter for importing cells and boundaries"
        
        if use_json:
            prereq = [args.xenium_ranger_assets]
        else:
            assert args.xenium_ranger_dir and os.path.exists(args.xenium_ranger_dir), "Please specify the --xenium-ranger-dir to import cells and boundaries"
            prereq = [os.path.join(args.xenium_ranger_dir, x) for x in [args.csv_cells, args.csv_boundaries, args.csv_clust, args.csv_diffexp] if x is not None]

        print(f" * Input JSON: {args.xenium_ranger_assets}" if use_json else f"Input Directory: {args.xenium_ranger_dir}")        
        import_cell_cmd=" ".join([
            "cartloader", "import_xenium_output",
            # actions
            "--all",
            f"--update-catalog" if not args.run_cartload2 and os.path.exists(catalog_yml) else "",
            # output
            f"--outprefix {cart_dir}/{args.cell_id}",
            # (conditional) input
            f"--in-json {args.xenium_ranger_assets}" if use_json else "",
            f"--in-dir {args.xenium_ranger_dir}" if not use_json else "",
            f"--csv-cells {args.csv_cells}" if (not use_json) and args.csv_cells else "",
            f"--csv-boundaries {args.csv_boundaries}" if (not use_json) and args.csv_boundaries else "",
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
            assert args.ome_tifs, "Please specify the --ome-tifs parameter to import images or set --detect-auto to automatically detect the images"

        assert args.image_ids, "Please specify at least of image ID --image-ids parameter to import images"
        if not args.image_colors:
            args.image_colors = colors[:len(args.image_ids)]
        assert len(args.image_ids) == len(args.image_colors), "Please specify the same number of image IDs and image colors in --image-ids and --image-colors parameters"
        
        for i in range(len(args.image_ids)):
            image_id = args.image_ids[i]
            image_color = args.image_colors[i].replace("#","")

            if use_json:
                prereq = [args.xenium_ranger_assets]
            else:
                prereq = [args.ome_tifs[i]]
                assert os.path.exists(args.ome_tifs[i]), f"Input OME TIFF file {args.ome_tifs[i]} does not exist. Please check the path."
            
            print(f" * Image ID: {image_id}", flush=True)
            print(f" * Image color: {image_color}", flush=True)

            import_image_cmd = " ".join([
                "cartloader", "import_image",
                # actions
                "--ome2png",
                "--png2pmtiles",
                f"--update-catalog" if not args.run_cartload2 and os.path.exists(catalog_yml) else "",
                # input
                f"--in-json {args.xenium_ranger_assets}" if use_json else "",
                f"--in-img {args.ome_tifs[i]}" if not use_json else "",
                # output
                f"--out-dir {cart_dir}",
                # id and color
                f"--img-id {image_id}",
                f"--colorize \"{image_color}\"",
                # aux params
                f"--upper-thres-quantile 0.95",
                f"--lower-thres-quantile 0.5",
                f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",

            ])
            import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["pmtiles", "gdaladdo"])
            import_image_cmd = add_param_to_cmd(import_image_cmd, args, ["restart","n_jobs"])
            
            run_command_w_preq(import_image_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)

            background_assets.append(os.path.join(cart_dir, f"{image_id}_assets.json"))
            background_pmtiles.append(os.path.join(cart_dir, f"{image_id}.pmtiles"))

    ficture_flags = []
    if args.run_ficture2:
        print("="*10, flush=True)
        print("Executing --run-ficture2 (execute via make)", flush=True)
        print("="*10, flush=True)
        # Note: new runs may be added to existing runs
        # In make, a rule runs only if the target is missing or older than its prerequisites — changing the commands alone doesn’t trigger re-execution.
        # So, should be a flag file relevant to all requested runs
        assert args.width, "Please specify the --width parameter for FICTURE2 analysis"
        assert args.n_factor, "Please specify the --n-factor parameter for FICTURE2 analysis"
        for tw in args.width.split(","):
            for nf in args.n_factor.split(","):
                ar=6
                fw=tw
                dc_prefix=os.path.join(fic_dir, f"t{tw}_f{nf}_p{fw}_a{ar}")
                ficture_flags.append(f"{dc_prefix}.done") # decoding step
                ficture_flags.append(f"{dc_prefix}_summary.done") # decoding info and de
                ficture_flags.append(f"{dc_prefix}.png") # decoding visual
                ficture_flags.append(os.path.join(fic_dir, f"ficture.params.json")) # decoding visual

        # prerequisites
        if args.sge_convert:
            prereq = [sge_flag]
        else:
            prereq = [sge_assets]
        
        print(f" * train width {args.width}", flush=True)
        print(f" * n factors {args.n_factor}", flush=True)
        ficture2_cmd = " ".join([
            "cartloader", "run_ficture2",
            # actions
            "--main",
            f"--out-dir {fic_dir}",
            f"--in-json {sge_assets}",
            f"--width {args.width}",
            f"--n-factor {args.n_factor}",
            f"--colname-feature {args.colname_feature}" if args.colname_feature else "",
            f"--colname-count {args.colname_count}" if args.colname_count else "",
        ])
        ficture2_cmd = add_param_to_cmd(ficture2_cmd, args, ["spatula", "ficture2"])
        ficture2_cmd = add_param_to_cmd(ficture2_cmd, args, ["restart","n_jobs","threads"])
        run_command_w_preq(ficture2_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)

    if args.run_cartload2:   
        print("="*10, flush=True)
        print("Executing --run-cartload2 (execute via make; function: run_cartload2_generic)", flush=True)
        print("="*10, flush=True)
        assert args.id is not None, "When --run-cartload2 is enabled, define the dataset ID by --id."

        # prerequisites
        prereq = [sge_flag]            
        if args.run_ficture2:
            print(f" * Importing FICTURE results: {fic_dir}", flush=True)
            prereq.extend(ficture_flags)
        elif args.import_ext_ficture2:
            print(f" * Importing external FICTURE results from {args.ext_fic_dir}", flush=True)
            assert args.ext_fic_dir and os.path.exists(args.ext_fic_dir), "Please specify the --ext-fic-dir parameter to import external FICTURE2 assets when --import-ext-ficture2 is enabled"
            ext_fic_params = os.path.join(args.ext_fic_dir, "ficture.params.json")
            prereq.append(ext_fic_params) 
        if args.import_cells:
            prereq.append(cell_assets)
        if args.import_images:
            prereq.extend(background_assets)

        #cmds = cmd_separator(cmds, f"Run cartload2 to create pmtiles for CartoScope")
        cartload_cmd = " ".join([
            "cartloader", "run_cartload2_generic",
            f"--out-dir {cart_dir}",
            f"--fic-dir {fic_dir}" if args.run_ficture2 else "",
            f"--ext-fic-dir {args.ext_fic_dir}" if args.import_ext_ficture2 else "",
            f"--sge_dir {args.sge_dir}" if (not args.run_ficture2) and (not args.import_ext_ficture2) else "",
            f"--cell-assets {cell_assets}" if args.import_cells else "",
            f"--background-assets {' '.join(background_assets)}" if background_assets else "",
            f"--id {args.id}" if args.id else "",
            f"--title {args.title}" if args.title else "",
            f"--desc {args.desc}" if args.desc else "",
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
        ])
        cartload_cmd = add_param_to_cmd(cartload_cmd, args, ["pmtiles", "gdaladdo", "spatula", "tippecanoe"])
        cartload_cmd = add_param_to_cmd(cartload_cmd, args, ["restart", "n_jobs", "threads"])
        run_command_w_preq(cartload_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)

    if args.upload_aws:
        print("="*10, flush=True)
        print("Executing --upload-aws (execute via make)", flush=True)
        print("="*10, flush=True)
        assert args.s3_bucket, "Specify a valid AWS S3 Bucket Name via --s3-bucket"

        prereq = [catalog_yml]

        aws_cmd = " ".join([
            "cartloader", "upload_aws",
            f"--in-dir {cart_dir}",
            f"--s3-dir s3://{args.s3_bucket}/{args.id}",
        ])
        aws_cmd = add_param_to_cmd(aws_cmd, args, ["aws"])
        aws_cmd = add_param_to_cmd(aws_cmd, args, ["restart", "n_jobs"])
        run_command_w_preq(aws_cmd, prerequisites=prereq, dry_run=args.dry_run, flush=True)

    zenodo_cartload_flag = f"cart_dir/cartload.zenodo.done"
    zenodo_basemap_flag = [f"{bg_pmtiles}.zenodo.done" for bg_pmtiles in background_pmtiles]
    if args.upload_zenodo:
        print("="*10, flush=True)
        print("Executing --upload-zenodo (execute directly)", flush=True)
        print("="*10, flush=True)
        assert args.zenodo_token and os.path.exists(args.zenodo_token), "Specify a valid file containing Zenodo token using --zenodo-token"
        if args.zenodo_title is None:
            args.zenodo_title = args.id
        assert args.zenodo_title, f"Specify a valid title for Zenodo deposition using --zenodo-title"
        prereq = [catalog_yml]

        zenodo_cmd = " ".join([
            "cartloader", "upload_zenodo",
            f"--in-dir {cart_dir}",
            f"--upload-method catalog",
            f"--zenodo-token {args.zenodo_token}",
            f"--zenodo-deposition-id {args.zenodo_deposition_id}" if args.zenodo_deposition_id else "",
            f"--title {args.zenodo_title}" if args.zenodo_title else "",
            f"--creators {args.creators}" if args.creators else "",
            f"--upload-type dataset" if not args.zenodo_deposition_id else "", # add type for new deposition
            f"--overwrite" if args.restart else ""
        ])
        if os.path.exists(zenodo_cartload_flag) and all(os.path.exists(f) for f in zenodo_basemap_flag) and not args.restart:
            print(f" * Skip --upload-zenodo since all flags exists. You can use --restart to force execute this step.", flush=True)
            print(import_cell_cmd, flush=True)
        else:
            run_command_w_preq(zenodo_cmd, prerequisites=[], dry_run=args.dry_run, flush=True)
    

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
