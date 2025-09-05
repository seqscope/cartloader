import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, re, inspect
import pandas as pd
import numpy as np

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, read_minmax, execute_makefile

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Convert a figure to pmtiles")

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run main commands (raster-dark, raster-light)')
    cmd_params.add_argument('--raster-dark', action='store_true', default=False, help='Create a dark-background raster PMTiles')
    cmd_params.add_argument('--raster-light', action='store_true', default=False, help='Create a light-background raster PMTiles')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Define the input file according to the user's needs.")
    inout_params.add_argument('--in-tsv', type=str, default=None, help='Input TSV file in micrometer scale')
    inout_params.add_argument('--colname-count', type=str, default='gn', help='Column name for feature counts')
    inout_params.add_argument('--in-minmax', type=str, default=None, help='Input minmax file in micrometer scale')
    inout_params.add_argument('--out-prefix', required=True, type=str, help='The output prefix. New directory will be created if needed')

    key_params = parser.add_argument_group("Key parameters")
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help='For the georeference and geo2tiff steps, define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--resample', type=str, default='cubic', help='For the geo2tiff step, define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    key_params.add_argument('--blocksize', type=int, default='512', help='For the geo2tiff step, define the blocksize (default: 512)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--spatula', type=str, default=f"{repo_dir}/submodules/spatula/bin/spatula", help='Path to spatula binary')
    aux_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    aux_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    aux_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate output files')
    aux_params.add_argument('--transparent-below', type=int, default=0, help='Threshold for transparent pixels below this value for dark background image (default: 0)')
    aux_params.add_argument('--transparent-above', type=int, default=255, help='Threshold for transparent pixels above this value for light background image (default: 255)')
    
    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default=None, help='The file name of the Makefile to generate (default: {out-prefix}.mk)')
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_tsv2mono(_args):

    args=parse_arguments(_args)

    if args.main:
        args.raster_dark = True
        args.raster_light = True

    # files
    out_dir = os.path.dirname(args.out_prefix)
    os.makedirs(out_dir, exist_ok=True)

    # start mm
    mm = minimake()

    ## check the existence of the applications
    scheck_app(args.spatula)
    scheck_app(args.pmtiles)
    scheck_app(args.gdal_translate)
    scheck_app(args.gdaladdo)
    # scheck_app(args.magick)

    ## check input tsv file header and get the column indices
    icol_x, icol_y, icol_cnt = None, None, None
    with gzip.open(args.in_tsv, 'rt') as f:
        header = f.readline().strip().split("\t")
        col2idx = {x:i for i,x in enumerate(header)}
        icol_x = col2idx["X"]
        icol_y = col2idx["Y"]
        if args.colname_count in col2idx:
            icol_cnt = col2idx[args.colname_count]

    ## get the minmax values
    minmax = read_minmax(args.in_minmax, "row")

    cmd_drawxy = " ".join([
        f"'{args.spatula}'",
        "draw-xy",
        "--tsv", args.in_tsv,
        "--icol-x", str(icol_x),
        "--icol-y", str(icol_y),
        "--icol-cnt", str(icol_cnt) if icol_cnt is not None else "-1",
        "--ullr", f"{minmax['xmin']},{minmax['ymin']},{minmax['xmax']},{minmax['ymax']}",
        "--auto-adjust",
        "--skip-lines", "1",
    ])

    cmds_dark  = cmd_separator([], f"Creating dark-background raster PMTiles: {args.out_prefix}-dark.pmtiles")
    cmds_light = cmd_separator([], f"Creating light-background raster PMTiles: {args.out_prefix}-light.pmtiles")

    ## create PNG files
    if args.transparent_below > 0:
        cmds_dark.append(f"{cmd_drawxy} --out {args.out_prefix}-dark-opaque.png")
        cmds_dark.append(" ".join([
            f"'{args.spatula}'",
            "png-mono2rgba",
            "--in", f"{args.out_prefix}-dark-opaque.png",
            "--out", f"{args.out_prefix}-dark.png",
            f"--transparent-below {args.transparent_below}" if args.transparent_below > 0 else ""
        ]))
    else:
        cmds_dark.append(f"{cmd_drawxy} --out {args.out_prefix}-dark.png")

    if args.transparent_above < 255:
        cmds_light.append(f"{cmd_drawxy} --invert --out {args.out_prefix}-light-opaque.png")
        cmds_light.append(" ".join([
            f"'{args.spatula}'",
            "png-mono2rgba",
            "--in", f"{args.out_prefix}-light-opaque.png",
            "--out", f"{args.out_prefix}-light.png",
            f"--transparent-above {args.transparent_above}" if args.transparent_above < 255 else ""
        ]))
    else:
        cmds_light.append(f"{cmd_drawxy} --invert --out {args.out_prefix}-light.png")

    ## create TIF files
    cmd_tif = " ".join([
        f"'{args.gdal_translate}'", 
        "-of", "GTiff", "-a_srs", args.srs,
        f"-a_ullr {minmax['xmin']} {minmax['ymin']} {minmax['xmax']} {minmax['ymax']}"
    ])
    cmds_dark.append(f"{cmd_tif} {args.out_prefix}-dark.png {args.out_prefix}-dark.pmtiles.tif")
    cmds_light.append(f"{cmd_tif} {args.out_prefix}-light.png {args.out_prefix}-light.pmtiles.tif")

    ## create MBTiles files
    cmd_mbt = " ".join([
        f"'{args.gdal_translate}'", 
        "-b 1 -b 2 -b 3 -b 4" if args.transparent_below > 0 else "-b 1",
        "-strict",
        "-co", "\"ZOOM_LEVEL_STRATEGY=UPPER\"",
        "-co", f"\"RESAMPLING={args.resample.upper()}\"",
        "-co", f"\"BLOCKSIZE={args.blocksize}\"",
        "-ot", "Byte", "-scale", "-of", "mbtiles",
        "-a_srs", args.srs
    ])
    cmds_dark.append(f"{cmd_mbt} {args.out_prefix}-dark.pmtiles.tif {args.out_prefix}-dark.pmtiles.mbtiles")
    cmd_mbt = " ".join([
        f"'{args.gdal_translate}'", 
        "-b 1 -b 2 -b 3 -b 4" if args.transparent_above < 255 else "-b 1",
        "-strict",
        "-co", "\"ZOOM_LEVEL_STRATEGY=UPPER\"",
        "-co", f"\"RESAMPLING={args.resample.upper()}\"",
        "-co", f"\"BLOCKSIZE={args.blocksize}\"",
        "-ot", "Byte", "-scale", "-of", "mbtiles",
        "-a_srs", args.srs
    ])
    cmds_light.append(f"{cmd_mbt} {args.out_prefix}-light.pmtiles.tif {args.out_prefix}-light.pmtiles.mbtiles")

    cmds_dark.append(f"'{args.gdaladdo}' {args.out_prefix}-dark.pmtiles.mbtiles -r {args.resample.lower()} 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536")
    cmds_light.append(f"'{args.gdaladdo}' {args.out_prefix}-light.pmtiles.mbtiles -r {args.resample.lower()} 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536")

    cmds_dark.append(f"'{args.pmtiles}' convert --force {args.out_prefix}-dark.pmtiles.mbtiles {args.out_prefix}-dark.pmtiles")
    cmds_light.append(f"'{args.pmtiles}' convert --force {args.out_prefix}-light.pmtiles.mbtiles {args.out_prefix}-light.pmtiles")

    if not args.keep_intermediate_files:
        cmds_dark.append(f"rm -f {args.out_prefix}-dark.pmtiles.mbtiles {args.out_prefix}-dark.pmtiles.tif {args.out_prefix}-dark.png {args.out_prefix}-dark-opaque.png")
        cmds_light.append(f"rm -f {args.out_prefix}-light.pmtiles.mbtiles {args.out_prefix}-light.pmtiles.tif {args.out_prefix}-light.png {args.out_prefix}-light-opaque.png")

    #cmds_dark.append(f"rm {args.out_prefix}-dark.pmtiles.mbtiles {args.out_prefix}-dark.pmtiles.tif")
    #cmds_light.append(f"rm {args.out_prefix}-light.pmtiles.mbtiles {args.out_prefix}-light.pmtiles.tif")

    cmds_dark.append(f"touch {args.out_prefix}-dark.pmtiles.done")
    cmds_light.append(f"touch {args.out_prefix}-light.pmtiles.done")

    if ( args.raster_dark ):
        mm.add_target(f"{args.out_prefix}-dark.pmtiles.done", [args.in_tsv, args.in_minmax], cmds_dark)
    if ( args.raster_light ):
        mm.add_target(f"{args.out_prefix}-light.pmtiles.done", [args.in_tsv, args.in_minmax], cmds_light)
    
    ## write makefile
    make_f=os.path.join(out_dir, args.makefn) if args.makefn is not None else f"{args.out_prefix}.mk"
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])