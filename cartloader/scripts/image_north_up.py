import sys, os, argparse, logging,  inspect, json, subprocess
import pandas as pd
import shutil
from cartloader.utils.utils import read_minmax
from cartloader.utils.image_helper import check_north_up

# get the path of the cu

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""
                                    Create a north-up tif from the input SGE file. 
                                    Input: can be provided by 1) provide input png file with minmax file, or 2) provide input tif file that is already georeferenced
                                     """)
    parser.add_argument('-i', '--in-png', type=str, help='Input png file')
    parser.add_argument('-t', '--in-tif', type=str, help='Input tif file')
    parser.add_argument('-o', '--out-tif', type=str, required=True, help='Output north-up tif file')
    parser.add_argument('--minmax', type=str, default=None, help='Input minmax file')
    parser.add_argument('--srs', type=str, default='EPSG:3857', help='Define the spatial reference system (default: EPSG:3857)')
    parser.add_argument('--resample', type=str, default='cubic', help='Define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    parser.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    parser.add_argument('--gdalwarp', type=str, default=f"gdalwarp", help='Path to gdalwarp binary')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Overwrite the output file if it exists')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def image_north_up(_args):
    args = parse_arguments(_args)

    assert args.in_png or args.in_tif, "Path not provided: --in-png or --in-tif"

    out_dir = os.path.dirname(args.out_tif)
    os.makedirs(out_dir, exist_ok=True)

    start_from_png=True
    if args.overwrite and args.in_png:
        print(f"--overwrite is enabled.")
        start_from_png=True
    if args.in_tif and os.path.exists(args.in_tif):
        start_from_png=False
    elif args.in_png:
        start_from_png=True
    else:
        raise ValueError(f"Please provide valid input files via --in-png ({args.in_png}) or --in-tif ({args.in_tif})")

    if not start_from_png:
        print(f"Starting with georeferenced tif: {args.in_tif}")
    else:
        print(f"Starting with png: {args.in_png}")
        # prepare in/out
        assert os.path.exists(args.in_png), f"File not found: {args.in_png} (--in-png)"
        assert args.minmax is not None, "Path not provided: --minmax"
        assert os.path.exists(args.minmax), f"File not found: {args.minmax} (--minmax)"
        if args.in_tif is None:
            args.in_tif = args.out_tif.replace(".tif", ".raw.tif")
        # convert png to georeferenced tif
        minmax=read_minmax(args.minmax, "row")
        cmd=f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr {minmax['xmin']} {minmax['ymin']} {minmax['xmax']} {minmax['ymax']} {args.in_png} {args.in_tif}"
        os.system(cmd)

    if check_north_up(args.in_tif):
        print(f"Image {args.in_tif} is already north-up.")
        # copy in_tif to args.out_tif
        shutil.copy(args.in_tif, args.out_tif)
    else:
        print(f"Image {args.in_tif} is not north-up. North-up the image.")
        # north-up the image
        cmd=f"{args.gdalwarp} -r {args.resample} -t_srs {args.srs} -overwrite -of GTiff {args.in_tif} {args.out_tif}"
        os.system(cmd)
    
    print(f"Image {args.out_tif} is now north-up.")


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
