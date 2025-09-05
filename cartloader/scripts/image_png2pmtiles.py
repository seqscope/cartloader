import sys, os, gzip, argparse, subprocess, inspect
import pandas as pd
import numpy as np

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, add_param_to_cmd, scheck_actions, write_dict_to_file, execute_makefile
from cartloader.utils.image_helper import orient2axisorder, update_orient

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
gdal_get_size_script = os.path.join(repo_dir, 'cartloader', "utils", "gdal_get_size.sh")

def cmds_for_dimensions(geotif_f, dim_f, gdalinfo="gdalinfo"):
    cmds = cmd_separator([], f"Extract dimensions from: {geotif_f}")
    dim_f = os.path.splitext(geotif_f)[0] + ".dim.tsv"
    cmds.append(f"{gdal_get_size_script} {geotif_f} {dim_f} {gdalinfo}")
    return cmds

def get_orientation_suffix(rotate, flip_vertical, flip_horizontal):
    ort_suffix = []
    if rotate is not None:
        ort_suffix.append(f"rot{rotate}")
    if flip_vertical:
        ort_suffix.append("vflip")
    if flip_horizontal:
        ort_suffix.append("hflip")
    return "-".join(ort_suffix)

def cmds_for_orientation(georef_f, dim_f,  ort_f, rotate, flip_vertical, flip_horizontal, mono=False, rgba=False):
    # extract axis order
    axis_order = orient2axisorder.get(
        (rotate, flip_vertical, flip_horizontal)
    )
    if axis_order is None:
        raise ValueError("Invalid combination of rotation and flip options.")

    # out dimensions
    if axis_order.startswith("1") or axis_order.startswith("-1"):
        out_dim="$WIDTH $HEIGHT"
    else:
        out_dim="$HEIGHT $WIDTH"
    
    # cmds
    msg=" ".join([f"Orientate the geotif, including ",
            "vertical flip;" if flip_vertical else "",
            "horizontal flip;" if flip_horizontal else "",
            f"rotate {rotate} degree clockwise;" if rotate is not None else "",
            georef_f])
    cmds = cmd_separator([], msg)
    cmds.append(f"WIDTH=$(awk '/WIDTH/' {dim_f}|cut -f 2) && \\")
    cmds.append(f"HEIGHT=$(awk '/HEIGHT/' {dim_f}|cut -f 2) && \\")
    cmd = " ".join([
            "gdalwarp",
            f'"{georef_f}"',  # Add quotes around file names to handle spaces
            f'"{ort_f}"',
            "-b 1" if mono else "-b 1 -b 2 -b 3 -b 4" if rgba else "-b 1 -b 2 -b 3",
            "-ct", f"\"+proj=pipeline +step +proj=axisswap +order={axis_order}\"",
            "-overwrite",
            "-ts", f"{out_dim}"
        ])
    cmds.append(cmd)
    return cmds
    
def create_mbtile_flag(mbtile_flag, mbtile_f, partial_db, journal_db):
    os.makedirs(os.path.dirname(mbtile_flag), exist_ok=True)
    conditions_met = (
        not os.path.exists(mbtile_flag) and
        os.path.exists(mbtile_f) and
        not os.path.exists(partial_db) and
        not os.path.exists(journal_db)
    )
    if conditions_met:
        try:
            with open(mbtile_flag, 'a'):
                pass 
            print(f"Created flag file: {mbtile_flag}")
        except OSError as e:
            print(f"Error creating flag file {mbtile_flag}: {e}")
            sys.exit(1)

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
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help='For --georeference and --geotif2mbtiles, define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--mono', action='store_true', default=False, help='Define if the input image is black-and-white and single-banded. Omit this if using --color-mode-record')
    key_params.add_argument('--rgba', action='store_true', default=False, help='RGBA, 4-banded image. Omit this if using --color-mode-record')
    key_params.add_argument('--color-mode-record', type=str, default=None, help='This argument is specifically designed to be used in "cartloader import_image"')
    key_params.add_argument('--resample', type=str, default='cubic', help='Resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    key_params.add_argument('--blocksize', type=int, default='512', help='Blocksize when creating mbtiles (default: 512)')
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
    # - mbtiles
    mbtile_f = f"{args.out_prefix}.pmtiles.mbtiles"
    mbtile_flag = f"{mbtile_f}.done"

    mbtile_f_ann=f"{args.out_prefix}.pmtiles.{args.resample}.mbtiles"

    # - pmtiles
    pmtiles_f=f"{args.out_prefix}.pmtiles"

    # args
    # - color arg (based on args.color_mode_record)
    if (args.flip_vertical or args.flip_horizontal or args.rotate is not None) or args.geotif2mbtiles:
        if args.color_mode_record:
            print(f"--color-mode-record is provided. Read color mode from {args.color_mode_record}")
            assert not args.rgba and not args.mono, "--color-mode-record cannot be used together with --mono or --rgba."
            assert os.path.exists(args.color_mode_record), f"File not found: {args.color_mode_record} (--color-mode-record)"
            with open(args.color_mode_record, 'r') as f:
                line = f.readline()
                color_mode = line.strip().lower()
                if color_mode not in {"rgb", "rgba", "mono"}:
                    raise ValueError(f"Invalid color mode '{color_mode}' in file: {args.color_mode_record}")
                if color_mode == "rgba":
                    args.rgba = True
                    args.mono = False
                elif color_mode == "mono":
                    args.rgba = False
                    args.mono = True
                elif color_mode == "rgb":
                    args.rgba = False
                    args.mono = False    

    # start mm
    mm = minimake()
    
    # 1. Georeferencing
    if args.georeference:
        scheck_app(args.gdal_translate)

        georef_f =f"{args.out_prefix}.georef.tif"
        cmds = cmd_separator([], f"Geo-referencing {args.in_img} to {georef_f}")

        if  args.georef_bounds is not None:
            ulx,uly,lrx,lry = args.georef_bounds.split(",")
        elif args.georef_bounds_tsv is not None:
            with open(args.georef_bounds_tsv, 'r') as f:
                line = f.readline()
            ullr_info = line.strip().lower()
            ulx,uly,lrx,lry = ullr_info.split(",")
        elif args.georef_pixel_tsv is not None:
            with gzip.open(args.georef_pixel_tsv, 'rt') as f:
                for i in range(3):
                    line = f.readline()
            ann2val = {x.split("=")[0]:x.split("=")[1] for x in line.strip().replace("##", "").split(";")}
            ulx = float(ann2val["OFFSET_X"]) # Upper left X-coordinate
            uly = float(ann2val["OFFSET_Y"]) # Upper left Y-coordinate
            lrx = float(ann2val["SIZE_X"]) + 1 + ulx # Lower Right X-coordinate
            lry = float(ann2val["SIZE_Y"]) + 1 + uly # Lower Right Y-coordinate 
        else:
            raise ValueError("Please provide bounds to georeference the image by --georef-bounds, --georef-bounds-tsv, --georef-pixel-tsv.")

        cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr {ulx} {uly} {lrx} {lry} {args.in_img} {georef_f}")
        mm.add_target(georef_f, [args.in_img], cmds)
    else:
        georef_f = args.in_img

    # 2. Orientation
    #  Check the orientation, and update it to the equivalent transformations
    if args.flip_vertical or args.flip_horizontal or args.rotate is not None:
        scheck_app(args.gdalinfo)
        
        args.rotate, args.flip_vertical, args.flip_horizontal = update_orient(args.rotate, args.flip_vertical, args.flip_horizontal)
        # Extract dimensions
        dim_f=georef_f.replace(".tif","") + ".dim.tsv"
        cmds=cmds_for_dimensions(georef_f, dim_f, args.gdalinfo)
        mm.add_target(dim_f, [georef_f], cmds)
    
        # Output FILENAME
        ort_suffix=get_orientation_suffix(args.rotate, args.flip_vertical, args.flip_horizontal)
        ort_f = f"{args.out_prefix}.{ort_suffix}.tif"

        # Orientation
        cmds = cmds_for_orientation(georef_f, dim_f, ort_f, args.rotate, args.flip_vertical, args.flip_horizontal, args.mono, args.rgba)
        mm.add_target(ort_f, [georef_f, dim_f], cmds)
    else:
        ort_f = georef_f

    # 3. Convert a raster image (GeoTIFF) into map tiles (MBTiles format)
    # * MBTiles expects top-left origin. If the input geotif use lower left was origin, GDAL will flip the origin to top left. 
    # * The origin could be find by running `gdalinfo <filename>` and check the `Pixel Size=(<x>,<y>)` field. If Pixel Y is negative, image is stored top-down. 
    # * Use a flag file instead of mbtiles_f to indicate the completion of the mbtiles conversion - to avoid silently failed or interrupted conversions.
    if args.geotif2mbtiles:
        scheck_app(args.gdal_translate)

        cmds = cmd_separator([], f"Converting from geotif to mbtiles: {ort_f}")

        ## Add cleanup step to prevent the case that the previous run was interrupted and left behind a .mbtiles file with a journal file or a partial_tiles.db file. Such unfinished jobs will not be availble to detect by the makefile.
        partial_db = mbtile_f.replace('.mbtiles', '.partial_tiles.db')
        journal_db = f"{mbtile_f}-journal"
        create_mbtile_flag(mbtile_flag, mbtile_f, partial_db, journal_db)
        cleanup_cmd = f"if [ -f {journal_db} ] || [ -f {partial_db} ] ; then echo 'Warning: Cleaning up incomplete previous conversion...' ; rm -f {mbtile_f} {journal_db} {partial_db} ; fi"
        cmds.append(cleanup_cmd)

        cmd = " ".join([
            args.gdal_translate, 
            "-b 1" if args.mono else "-b 1 -b 2 -b 3 -b 4" if args.rgba else "-b 1 -b 2 -b 3", # Use the Red, Green, Blue bands
            "-strict",                                      # Enforce strict format compliance
            "-co", "\"ZOOM_LEVEL_STRATEGY=UPPER\"",         # Use higher zoom levels to preserve detail
            "-co", f"\"RESAMPLING={args.resample}\"",       # Use cubic interpolation when scaling tiles
            "-co", f"\"BLOCKSIZE={args.blocksize}\"",       # Use 512x512 pixel tile blocks (standard for MBTiles)
            "-ot", "Byte",                                  # Convert pixel values to 8-bit integers (0-255)
            "-scale" if args.mono else "",                  # Automatically scale pixel values to 0-255
            "-of", "mbtiles",
            "-a_srs", args.srs,                             # default to EPSG:3857. Assign the target projection: Web Mercator (used by web maps)
            ort_f, 
            mbtile_f
        ])
        cmds.append(cmd)
        
        validation_cmd = f" [ -f {mbtile_f} ]  && [ ! -f {journal_db} ] && [ ! -f {partial_db} ] && touch {mbtile_flag}"
        cmds.append(validation_cmd)

        mm.add_target(mbtile_flag, [ort_f], cmds)
    
    # 4. Convert mbtiles to pmtiles
    if args.mbtiles2pmtiles:
        scheck_app(args.pmtiles)
        scheck_app(args.gdaladdo)

        cmds = cmd_separator([], f"Resampling mbtiles and converting to pmtiles: {mbtile_f}")
        cmds.append(f"cp {mbtile_f} {mbtile_f_ann}")
        cmds.append(f"'{args.gdaladdo}' {mbtile_f_ann} -r {args.resample} 2 4 8 16 32 64 128 256") # Build internal overviews. 
        cmds.append(f"'{args.pmtiles}' convert --force {mbtile_f_ann} {pmtiles_f}")
        cmds.append(f" [ -f {pmtiles_f} ] && rm {mbtile_f_ann}") # clean temp files
        mm.add_target(pmtiles_f, [mbtile_flag], cmds)

    ## write makefile
    make_f=os.path.join(out_dir, args.makefn) if args.makefn is not None else f"{args.out_prefix}.png2pmtiles.mk"
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)

    ## write a json file 
    img_id = os.path.basename(args.out_prefix)
    image_assets={
        img_id: pmtiles_f
    }
    write_dict_to_file(image_assets, f"{args.out_prefix}_assets.json", check_equal=True)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
