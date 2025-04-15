import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess
import rasterio # for extracting bounds

from cartloader.utils.utils import cmd_separator
from cartloader.utils.minimake import minimake
from cartloader.utils.image_helper import update_orient

from cartloader.scripts.run_fig2pmtiles import get_orientation_suffix, cmds_for_dimensions, cmds_for_orientation

# get the current path
current_path = os.path.realpath(__file__)
cartloader_dir=os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
gdal_get_size_script = os.path.join(cartloader_dir, 'cartloader', "utils", "gdal_get_size.sh")

def hist_stitch(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="""
        Stitching multiple tif tiles into one tif. 
        If georeferencing is required, there are two options: 1) set georef=True with ullr in local (per-tile) coordinates and provide --in-offsets; 2) set georef=True with ullr in global coordinates and set --in-offsets=None.
        If the tiles are already georeferenced with local coordinates, use --in-offsets to adjust positions.
        If the tiles are already georeferenced with global coordinates, set georef=False and --in-offsets=None.
        """
    )
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Simulate the process without executing commands (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore all intermediate files and start from the beginning (default: False)')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default=None, help='Makefile name. By default, it will be named as hist_stitch_<filename>.mk based on the output file name.')
    run_params.add_argument('--threads', type=int, default=1, help='Maximum number of threads to use in each process (default: 1)')
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/Output Parameters")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of input tiles in the format: <path>,<row>,<col>,<georef>,<bounds>,<rotate>,<vertical_flip>,<horizontal_flip>. "
                                                                                "<georef> is a boolean indicating if georeferencing is needed. "
                                                                                "<bounds> is '<ulx>,<uly>,<lrx>,<lry>' if georeferencing is required. "
                                                                                "<rotate> is one of 0, 90, 180, 270. <vertical_flip> and <horizontal_flip> are booleans indicating if the image should be flipped vertically or horizontally.")
    inout_params.add_argument("--output", type=str, help="Output path for the stitched tif file.")
    inout_params.add_argument("--in-offsets", type=str, default=None, help="Path to the input offsets file, in which the following columns are required: row, col, x_offset, y_offset, units_per_um")
    
    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    env_params.add_argument('--gdalbuildvrt', type=str, default=f"gdalbuildvrt", help='Path to gdalbuildvrt binary')
    env_params.add_argument('--srs', type=str, default='EPSG:3857', help='Spatial reference system (default: EPSG:0)')
    args = parser.parse_args(_args)

    mm = minimake()

    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)

    out_bn = os.path.basename(args.output).replace(".tif", "")
    if args.makefn is None:
        args.makefn = f"hist_stitch_{out_bn}.mk"

    # read in_offsets
    idx2offsets = {}
    if args.in_offsets is not None:
        offsets = pd.read_csv(args.in_offsets, sep="\t")

        # assure the type
        offsets["row"] = offsets["row"].astype(int).astype(str)
        offsets["col"] = offsets["col"].astype(int).astype(str)
        offsets["x_offset"] = offsets["x_offset"].astype(float)
        offsets["y_offset"] = offsets["y_offset"].astype(float)
        offsets["units_per_um"] = offsets["units_per_um"].astype(float)

        # assuming the coordinates in hist in um
        offsets["x_offset_um"] = offsets["x_offset"] / offsets["units_per_um"]
        offsets["y_offset_um"] = offsets["y_offset"] / offsets["units_per_um"]
   
        # create dict between (row, col) to x_offset, y_offset
        idx2offsets = {
            (row["row"], row["col"]): (row["x_offset_um"], row["y_offset_um"])
            for _, row in offsets.iterrows()
        }

    tiles = []
    for in_tile in args.in_tiles:
        # input params
        tif, row, col, georef, georef_bound, rotate, vflip, hflip = in_tile.split(",")

        assert os.path.exists(tif), f"Input histology file {tif} (index: {row}, {col}) does not exist."
        row, col = str(int(row)), str(int(col))
        georef = str(georef).lower() == "true"
        georef_bound = georef_bound if georef_bound else None
        rotate = str(rotate) if rotate else None
        vflip = str(vflip).lower() == "true"
        hflip = str(hflip).lower() == "true"
        
        # tile_prefix for orientation
        tile_prefix = os.path.join(out_dir, os.path.basename(tif).replace(".tif", ""))

        # Update to global coordinates
        # 1) georef 
        x_offset, y_offset = idx2offsets.get((row, col), (0, 0))
        
        # - tif without coordinates with/wo offsets (update to global coordinates)
        if georef: 
            assert georef_bound is not None, "Georeferencing requires bounds. Format: <ulx>,<uly>,<lrx>,<lry>"
            ulx, uly, lrx, lry = map(float, georef_bound.split(","))
            ulx = ulx + x_offset 
            uly = uly + y_offset
            lrx = lrx + x_offset
            lry = lry + y_offset
            ullr = f"{ulx} {uly} {lrx} {lry}"
        # - tif with local coordinates (add offsets to update to global coordinates)
        elif x_offset != 0 or y_offset != 0:
            # extract bounds from the tif and add offsets from gdalinfo
            with rasterio.open(tif) as src:
                bounds = src.bounds             # returns: BoundingBox(left=min_x, bottom=min_y, right=max_x, top=max_y)
            ulx = bounds.left + x_offset
            uly = bounds.top + y_offset
            lrx = bounds.right + x_offset
            lry = bounds.bottom + y_offset
            ullr = f"{ulx} {uly} {lrx} {lry}"
        # - tif with global coordinates (no offsets)
        else:
            ullr = None

        if ullr is not None:
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Georeferencing with ullr ({ullr}). Input params: georef={georef}, bounds={georef_bound}, offsets=({x_offset}, {y_offset})")
            georef_f=f"{tile_prefix}.georef.tif"
            cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr {ullr} {tif} {georef_f}")
            mm.add_target(georef_f, [tif], cmds)
        else:
            georef_f = tif
        
        # 2) rotate
        if rotate is not None or vflip or hflip:
            rotate, vflip, hflip = update_orient(rotate, vflip, hflip, tif)
            # Generate Dimension file to provide height and width
            dim_f = georef_f.replace(".tif", "") + ".dim.tsv"
            cmds = cmds_for_dimensions(georef_f, dim_f)
            mm.add_target(dim_f, [georef_f], cmds)
            # Generate the orientated file
            ort_suffix = get_orientation_suffix(rotate, vflip, hflip)
            ort_f = f"{tile_prefix}.{ort_suffix}.tif"
            cmds = cmds_for_orientation( georef_f, dim_f, ort_f, rotate, vflip, hflip)
            mm.add_target(ort_f, [georef_f, dim_f], cmds) 
        else:
            ort_f = georef_f

        tiles.append(ort_f)

    # 3) stitch
    # gdalbuildvrt + gdal_translate
    cmds= cmd_separator([], f"Stitching the input into {args.output}")
    out_vrt = args.output.replace(".tif", ".vrt")
    cmd = " ".join([args.gdalbuildvrt, out_vrt] + tiles)
    cmds.append(cmd)
    cmd = " ".join([
        args.gdal_translate,
        "-of GTiff",
        out_vrt,
        args.output])
    cmds.append(cmd)
    mm.add_target(args.output, tiles, cmds)
            
    # write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    
    make_f = os.path.join(out_dir, args.makefn)
    mm.write_makefile(make_f)
    if args.dry_run:
        dry_cmd=f"make -f {make_f} -n {'-B' if args.restart else ''} "
        os.system(dry_cmd)
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}")
    else:
        exe_cmd=f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}"
        result = subprocess.run(exe_cmd, shell=True)
        if result.returncode != 0:
            print(f"Error in executing: {exe_cmd}")
            sys.exit(1)

if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
