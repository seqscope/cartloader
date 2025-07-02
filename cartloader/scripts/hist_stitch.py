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
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/Output Parameters")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of input tiles in the format: <row>,<col>,<path>,<georef>,<georef_tsv>,<georef_bounds>,<rotate>,<vertical_flip>,<horizontal_flip>. "
                                                                                "<georef> is a boolean indicating if georeferencing is needed. "
                                                                                "<georef_tsv> is the path to the *.pixel.sorted.tsv.gz from run_ficture to provide georeferenced bounds."
                                                                                "<georef_bounds> is '<ulx>_<uly>_<lrx>_<lry>' if georeferencing is required. "
                                                                                "<rotate> is one of 0, 90, 180, 270. <vertical_flip> and <horizontal_flip> are booleans indicating if the image should be flipped vertically or horizontally.")
    inout_params.add_argument("--output", type=str, help="Output path for the stitched tif file.")
    inout_params.add_argument("--in-offsets", type=str, default=None, help="Path to the input offsets file, in which the following columns are required: row, col, x_offset, y_offset, units_per_um")
    inout_params.add_argument("--crop-tile-by-minmax", action='store_true', default=False, help="Crop the tile by minmax.")
    inout_params.add_argument("--in-minmax", type=str, default=None, help="Required if --crop-tile-by-minmax is enabled. Path to the input minmax file, in which the following columns are required: row, col, global_xmin_um, global_ymin_um, global_xmax_um, global_ymax_um")
    inout_params.add_argument("--downsize", action='store_true', default=False, help="Downsize the image to a proportion of the original size. This is only for testing purposes.")
    inout_params.add_argument("--downsize-prop", type=float, default=0.25, help="Downsize proportion (default: 0.25).")
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
        offsets["x_offset_unit"] = offsets["x_offset_unit"].astype(float)
        offsets["x_offset_unit"] = offsets["x_offset_unit"].astype(float)
        offsets["units_per_um"] = offsets["units_per_um"].astype(float)

        # assuming the coordinates in hist in um
        offsets["x_offset_um"] = offsets["x_offset_unit"] / offsets["units_per_um"]
        offsets["y_offset_um"] = offsets["y_offset_unit"] / offsets["units_per_um"]
   
        # create dict between (row, col) to x_offset, y_offset
        idx2offsets = {
            (row["row"], row["col"]): (row["x_offset_um"], row["y_offset_um"])
            for _, row in offsets.iterrows()
        }
    if args.crop_tile_by_minmax:
        minmax = pd.read_csv(args.in_minmax, sep="\t")
        # assure the type
        minmax["row"] = minmax["row"].astype(int).astype(str)
        minmax["col"] = minmax["col"].astype(int).astype(str)
        minmax["global_xmin_um"] = minmax["global_xmin_um"].astype(float)
        minmax["global_ymin_um"] = minmax["global_ymin_um"].astype(float)
        minmax["global_xmax_um"] = minmax["global_xmax_um"].astype(float)
        minmax["global_ymax_um"] = minmax["global_ymax_um"].astype(float)
        # create dict between (row, col) to xmin, ymin, xmax, ymax
        idx2minmax = {
            (row["row"], row["col"]): (row["global_xmin_um"], row["global_ymin_um"], row["global_xmax_um"], row["global_ymax_um"])
            for _, row in minmax.iterrows()
        }
    
    tiles = []
    for in_tile in args.in_tiles:
        # input params
        row, col, tif, georef, georef_tsv, georef_bound, rotate, vflip, hflip = in_tile.split(",")

        assert os.path.exists(tif), f"Input histology file {tif} (index: {row}, {col}) does not exist."
        tile_prefix = os.path.join(out_dir, os.path.basename(tif).replace(".tif", "")) # tile_prefix for intermediate files

        row, col = str(int(row)), str(int(col))

        # actions
        # - georef
        georef = str(georef).lower() == "true"
        georef_tsv = georef_tsv if georef_tsv else None
        georef_bound = georef_bound if georef_bound else None
        # - orient
        rotate = str(rotate) if rotate else None
        vflip = str(vflip).lower() == "true"
        hflip = str(hflip).lower() == "true"
        rotate, vflip, hflip = update_orient(rotate, vflip, hflip, tif)
        orient= True if rotate is not None or vflip or hflip else False
        # - add offsets
        x_offset, y_offset = idx2offsets.get((row, col), (0, 0))
        add_offsets = True if x_offset != 0 or y_offset != 0 else False

        current_f = tif
        
        # 1) georef (tif without coordinates)
        if georef: 
            assert not (georef_bound is None and georef_tsv is None), f"Tile (row {row}, col {col}): Georeferencing requires either georef_bound or georef_tsv to be provided."
            assert not (georef_bound is not None and georef_tsv is not None), f"Tile (row {row}, col {col}): Georeferencing requires only one of georef_bound or georef_tsv to be provided, not both."
            if georef_bound is not None:
                ulx, uly, lrx, lry = map(float, georef_bound.split("_"))
            elif georef_tsv is not None:
                with gzip.open(georef_tsv, 'rt') as f:
                    for i in range(3):
                        line = f.readline()
                ann2val = {x.split("=")[0]:x.split("=")[1] for x in line.strip().replace("##", "").split(";")}
                ulx = float(ann2val["OFFSET_X"]) # Upper left X-coordinate
                uly = float(ann2val["OFFSET_Y"]) # Upper left Y-coordinate
                lrx = float(ann2val["SIZE_X"])+1+ulx # Lower Right X-coordinate
                lry = float(ann2val["SIZE_Y"])+1+uly # Lower Right Y-coordinate 
            ullr=f"{ulx} {uly} {lrx} {lry}"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Georeferencing with ullr ({ullr})")
            georef_f=f"{tile_prefix}.georef.tif"
            cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr {ullr} {tif} {georef_f}")
            mm.add_target(georef_f, [tif], cmds)
            current_f = georef_f

        # 2) rotate
        if orient:
            # Generate Dimension file to provide height and width
            # dim_f = current_f.replace(".tif", "") + ".dim.tsv"
            dim_f = os.path.splitext(current_f)[0] + ".dim.tsv"     # Modify this since encountered a case with two tif in the figure name
            cmds = cmds_for_dimensions(current_f, dim_f)
            mm.add_target(dim_f, [current_f], cmds)
            # Generate the orientated file
            ort_suffix = get_orientation_suffix(rotate, vflip, hflip)
            ort_f = f"{tile_prefix}.{ort_suffix}.tif"
            cmds = cmds_for_orientation( current_f, dim_f, ort_f, rotate, vflip, hflip)
            mm.add_target(ort_f, [current_f, dim_f], cmds) 
            current_f = ort_f
            # update the ullr
            
        # 3) add offsets
        if add_offsets:
            # extract local bounds from the ort_f
            #ullr_f = current_f.replace(".tif", f".bounds.offsets_x{x_offset}_y{y_offset}.txt")
            ullr_f =f"{tile_prefix}.globalcoord.ullr.txt"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Extracting bounds from {current_f} and adding offsets ({x_offset}, {y_offset}).")
            cmds.append(f"cartloader image_extract_ullr --input {current_f} --output {ullr_f} --x_offset {x_offset} --y_offset {y_offset}")
            mm.add_target(ullr_f, [current_f], cmds)
            # read the bounds
            gcoord_f = f"{tile_prefix}.globalcoord.tif"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Georeferencing with ullr after adding offsets ({x_offset}, {y_offset}).")
            # append a bash command to get the ullr from the first row in ullr_f
            cmds.append(f"ULLR=$(head -n 1 {ullr_f}) && \\")
            cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr $ULLR {current_f} {gcoord_f}")
            mm.add_target(gcoord_f, [current_f, ullr_f], cmds)
            current_f = gcoord_f
        
        # 4) crop
        if args.crop_tile_by_minmax:
            # read the minmax
            xmin, ymin, xmax, ymax = idx2minmax.get((row, col), (0, 0, 0, 0))
            assert not (xmin == 0 and ymin == 0 and xmax == 0 and ymax == 0), f"Tile (row {row}, col {col}): The minmax for the tile is not provided in the input minmax file."
            ullr_f = f"{tile_prefix}.crop.ullr.txt"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Extracting bounds from {current_f} and cropping the image by minmax ({xmin}, {ymin}, {xmax}, {ymax}).")
            cmds.append(f"cartloader image_extract_ullr --input {current_f} --output {ullr_f} --crop-by-minmax --minmax {xmin},{xmax},{ymin},{ymax}")
            mm.add_target(ullr_f, [current_f], cmds)
            # crop the image
            crop_f = f"{tile_prefix}.crop.tif"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Cropping the image by minmax ({xmin}, {ymin}, {xmax}, {ymax}).")
            cmds.append(f"ULLR=$(head -n 1 {ullr_f}) && \\")
            cmds.append(f"{args.gdal_translate} -of GTiff -projwin $ULLR {current_f} {crop_f}")
            mm.add_target(crop_f, [current_f, ullr_f], cmds)
            current_f = crop_f

        tiles.append(current_f)

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
        
    # 4) downsize
    if args.downsize:
        downsize_f = args.output.replace(".tif", f".downsize_{args.downsize_prop}.tif")
        cmds = cmd_separator([], f"Downsizing the image to {args.downsize_prop} of the original size.")
        downsize_percentage = f"{args.downsize_prop * 100}%"
        cmds.append(f"gdal_translate -outsize {downsize_percentage} {downsize_percentage} -r cubic {args.output} {downsize_f}")
        mm.add_target(downsize_f, [args.output], cmds)

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
