import sys, os, argparse, inspect

from cartloader.image import configure_color_mode, register_png2pmtiles_pipeline
from cartloader.utils.minimake import minimake
from cartloader.utils.utils import scheck_actions, write_dict_to_file, execute_makefile
from cartloader.utils.orient_helper import update_orient

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="""
                                    Converts a TIF or PNG image to PMTiles by first generating MBTiles and then converting to PMTiles. 
                                    Supports optional georeferencing, rotation, and flipping of the input image.
                                     """)

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default=None, help='The file name of the Makefile to generate (default: {out-prefix}.png2pmtiles.mk)')

    cmd_params = parser.add_argument_group("Commands")
    cmd_params.add_argument('--geotif2mbtiles', action='store_true', default=False, help="Convert from a georeferenced file to mbtiles")
    cmd_params.add_argument('--mbtiles2pmtiles', action='store_true', default=False, help="Convert from mbtile to pmtile")
    cmd_params.add_argument('--georeference', action='store_true', default=False, help='Create a geotiff file from PNG or TIF file. If enabled, the user must provide georeferenced bounds using --in-tsv or --in-bounds.')
    cmd_params.add_argument('--rotate', type=str, default=None, choices=["90", "180", "270"],  help='Rotate the image by 90, 180, or 270 degrees clockwise. Rotate precedes flip.')
    cmd_params.add_argument('--flip-vertical', action='store_true', default=False, help='Flip the image vertically (flipped along Y-axis). Rotate precedes flip.')
    cmd_params.add_argument('--flip-horizontal', action='store_true', default=False, help='Flip the image horizontally (flipped along X-axis). Rotate precedes flip.')

    inout_params = parser.add_argument_group("Input/Output Parameters")
    inout_params.add_argument('--in-img', type=str, help='Path to the input image file (PNG or TIF)')
    inout_params.add_argument('--out-prefix', type=str, help='Prefix for naming the output files.')

    key_params = parser.add_argument_group("Key parameters")
    key_params.add_argument('--georef-pixel-tsv', type=str, default=None, help='If --georeference is required, use the *.pixel.sorted.tsv.gz from run_ficture to provide georeferenced bounds.')
    key_params.add_argument('--georef-bounds-tsv', type=str, default=None, help='If --georeference is required, provide the bounds via a tsv file. This TSV should include 1 line with <ulx>,<uly>,<lrx>,<lry> ')
    key_params.add_argument('--georef-bounds', type=str, default=None, help='If --georeference is required, provide the bounds in the format of "<ulx>,<uly>,<lrx>,<lry>", which represents upper-left X, upper-left Y, lower-right X, lower-right Y.')
    key_params.add_argument('--georef-detect', type=str, default=None, help='If --georeference is required, detect the bounds from the image metadata. Supported formats: ome')
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help='For --georeference and --geotif2mbtiles, define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--mono', action='store_true', default=False, help='Define if the input image is black-and-white and single-banded. Omit this if using --color-mode-record')
    key_params.add_argument('--rgba', action='store_true', default=False, help='RGBA, 4-banded image. Omit this if using --color-mode-record')
    key_params.add_argument('--color-mode-record', type=str, default=None, help='This argument is specifically designed to be used in "cartloader import_image"')
    key_params.add_argument('--resample', type=str, default='cubic', help='Resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    key_params.add_argument('--blocksize', type=int, default=512, help='Blocksize when creating mbtiles (default: 512)')
    #key_params.add_argument('--remove-intermediate-files', action='store_true', default=False, help='If set, remove intermediate files (e.g., .mbtiles) after generating the final output.')
    
    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles (default: pmtiles)')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binar (default: gdaladdo)')
    env_params.add_argument('--gdalinfo', type=str, default=f"gdalinfo", help='Path to gdalinfo binary (default: gdalinfo)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def image_png2pmtiles(_args):

    args=parse_arguments(_args)

    # actions
    scheck_actions(args, ["--geotif2mbtiles", "--mbtiles2pmtiles", "--georeference", "--rotate", "--flip-vertical", "--flip-horizontal"], context="actions")

    # dirs/files
    # - output dir
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir != "":
        os.makedirs(out_dir, exist_ok=True)
    
    # - input image
    assert args.in_img is not None, f"Path not provided: --in-img"
    assert os.path.exists(args.in_img), f"File not found: {args.in_img} (--in-img)"

    # args
    if args.color_mode_record:
        print(f"--color-mode-record is provided. Read color mode from {args.color_mode_record}")
    configure_color_mode(args)

    if args.flip_vertical or args.flip_horizontal or args.rotate is not None:
        args.rotate, args.flip_vertical, args.flip_horizontal = update_orient(
            args.rotate, args.flip_vertical, args.flip_horizontal
        )

    # start mm
    mm = minimake()
    result = register_png2pmtiles_pipeline(mm, args)

    ## write makefile
    make_f=os.path.join(out_dir, args.makefn) if args.makefn is not None else f"{args.out_prefix}.png2pmtiles.mk"
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)

    ## write a json file 
    img_id = os.path.basename(args.out_prefix)
    if result.pmtiles_path:
        image_assets={
            img_id: result.pmtiles_path
        }
        write_dict_to_file(image_assets, f"{args.out_prefix}_assets.json", check_equal=True)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
