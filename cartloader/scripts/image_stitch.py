import sys, os, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess
from types import SimpleNamespace

from cartloader.utils.utils import cmd_separator, execute_makefile
from cartloader.utils.minimake import minimake
from cartloader.utils.orient_helper import update_orient
from cartloader.image import register_georeference_stage, register_orientation_stage

# get the current path
current_path = os.path.realpath(__file__)
cartloader_dir=os.path.dirname(os.path.dirname(os.path.dirname(current_path)))


def _parse_georef_info(georef_info: str, row: str, col: str, args) -> SimpleNamespace:
    if not georef_info or not georef_info.strip():
        raise ValueError(
            f"Tile (row {row}, col {col}): georeferencing requested but georef_info is empty"
        )

    if ";" in georef_info:
        raise ValueError(
            f"Tile (row {row}, col {col}): invalid georef_info '{georef_info}'. Provide a single key:value entry"
        )

    token = georef_info.strip()
    if not token:
        raise ValueError(
            f"Tile (row {row}, col {col}): georeferencing requested but georef_info is empty"
        )

    allowed_keys = ("bounds", "bounds_tsv", "bounds_pixel_tsv", "detect")
    token_lower = token.lower()
    key = None
    value = None
    for candidate in allowed_keys:
        prefix_colon = f"{candidate}:"
        prefix_eq = f"{candidate}="
        if token_lower.startswith(prefix_colon):
            key = candidate
            value = token[len(prefix_colon):]
            break
        if token_lower.startswith(prefix_eq):
            key = candidate
            value = token[len(prefix_eq):]
            break

    if key is None:
        raise ValueError(
            f"Tile (row {row}, col {col}): invalid georef_info '{token}'. Expected one of {allowed_keys} followed by ':'"
        )

    value = value.strip()
    if not value:
        raise ValueError(
            f"Tile (row {row}, col {col}): georef_info '{token}' is missing a value"
        )

    geo_args = SimpleNamespace(
        gdal_translate=args.gdal_translate,
        srs=args.srs,
        georef_bounds=None,
        georef_bounds_tsv=None,
        georef_pixel_tsv=None,
        georef_detect=None,
    )

    if key == "bounds":
        geo_args.georef_bounds = value.replace("_", ",")
    elif key == "bounds_tsv":
        geo_args.georef_bounds_tsv = value
    elif key == "bounds_pixel_tsv":
        geo_args.georef_pixel_tsv = value
    elif key == "detect":
        geo_args.georef_detect = value

    return geo_args

def image_stitch(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="""
        Stitching multiple TIFF tiles into one TIFF file. 
        If georeferencing is required (set georef=True), choose one of the following options to provide global coordinates:
            1) Provide local (per-tile) coordinates as ullr using <georef_tsv> or <georef_bounds>, and use --in-offsets to convert them to global coordinates.
            2) Provide global coordinates directly as ullr using <georef_tsv> or <georef_bounds>, and omit --in-offsets.
        If the tiles are already georeferenced (georef is not needed; set georef=False):
            1.	If georeferenced using local coordinates, use --in-offsets to adjust positions.
            2.	If georeferenced using global coordinates, set georef=False and omit --in-offsets.
        """
    )
    run_params = parser.add_argument_group("Run Options", "Run options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Simulate the process without executing commands')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', '-j', type=int, default=1, help='Number of jobs (processes) to run in parallel (default: 1)')
    run_params.add_argument('--makefn', type=str, default=None, help='Makefile name. By default, it will be named as image_stitch_<filename>.mk based on the output file name.')
    
    inout_params = parser.add_argument_group("Input/Output Parameters", "Input/Output Parameters")
    inout_params.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of input tiles in the format: <row>,<col>,<path>,<georef>,<georef_info>,<rotate>,<vertical_flip>,<horizontal_flip>."
                                                                                "<georef>: Boolean flag indicating whether georeferencing is required."
                                                                                "<georef_info>: key:value entry describing bounds; supported keys: bounds (ulx_uly_lrx_lry), bounds_tsv, bounds_pixel_tsv, detect; supported detect type: ome."
                                                                                "<rotate>: Rotation angle (0, 90, 180, or 270)."
                                                                                "<vertical_flip> and <horizontal_flip> are booleans indicating if the image should be flipped vertically or horizontally.")
    inout_params.add_argument("--output", type=str, help="Path to the output stitched TIFF file.")
    inout_params.add_argument("--in-offsets", type=str, default=None, help="Path to the offsets file. Must contain columns: row, col, x_offset, y_offset, and units_per_um.")
    inout_params.add_argument("--crop-tile-by-minmax", action='store_true', default=False, help="Enable cropping tiles using min/max bounds (default: False).")
    inout_params.add_argument("--in-minmax", type=str, default=None, help="Required if --crop-tile-by-minmax is set. Path to the input min/max file with columns: row, col, global_xmin_um, global_ymin_um, global_xmax_um, global_ymax_um.")
    inout_params.add_argument("--downsize", action='store_true', default=False, help="Enable Downsizing the stitched image for testing purposes (default: False).")
    inout_params.add_argument("--downsize-prop", type=float, default=0.25, help="Proportion to downsize the image by (default: 0.25).")
    
    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdalbuildvrt', type=str, default=f"gdalbuildvrt", help='Path to gdalbuildvrt binary (default: gdalbuildvrt)')
    env_params.add_argument('--gdalinfo', type=str, default='gdalinfo', help='Path to gdalinfo binary (default: gdalinfo)')
    env_params.add_argument('--gdalwarp', type=str, default='gdalwarp', help='Path to gdalwarp binary (default: gdalwarp)')
    env_params.add_argument('--srs', type=str, default='EPSG:3857', help='Spatial reference system (default: EPSG:3857)')

    color_params = parser.add_argument_group("Color Options", "Band configuration for orientation")
    color_params.add_argument('--mono', action='store_true', default=False, help='Treat tiles as single-band images during orientation')
    color_params.add_argument('--rgba', action='store_true', default=False, help='Treat tiles as 4-band RGBA images during orientation (mutually exclusive with --mono)')
    args = parser.parse_args(_args)

    if args.mono and args.rgba:
        parser.error("--mono and --rgba are mutually exclusive")

    mm = minimake()

    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)

    out_bn = os.path.basename(args.output).replace(".tif", "")
    if args.makefn is None:
        args.makefn = f"image_stitch_{out_bn}.mk"

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

        # assuming the coordinates in the image in um
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
        parts = in_tile.split(",")
        if len(parts) != 8:
            raise ValueError(
                "--in-tiles entries must have 8 fields (row,col,path,georef,georef_info,rotate,vertical_flip,horizontal_flip)"
            )
        row, col, tif, georef, georef_info, rotate, vflip, hflip = parts

        assert os.path.exists(tif), f"Tile (row {row}, col {col}): file not found {tif} (--in-tiles)"
        tile_prefix = os.path.join(out_dir, os.path.basename(tif).replace(".tif", "")) # tile_prefix for intermediate files

        row, col = str(int(row)), str(int(col))

        # actions
        # - georef
        georef = str(georef).lower() == "true"
        georef_info = georef_info.strip()
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
            geo_args = _parse_georef_info(georef_info, row, col, args)

            current_f = register_georeference_stage(
                mm,
                geo_args,
                in_img=current_f,
                out_prefix=tile_prefix,
            )

        # 2) rotate
        if orient:
            orient_args = SimpleNamespace(
                rotate=rotate,
                flip_vertical=vflip,
                flip_horizontal=hflip,
                mono=args.mono,
                rgba=args.rgba,
                gdalinfo=args.gdalinfo,
                gdalwarp=args.gdalwarp,
            )

            current_f = register_orientation_stage(
                mm,
                orient_args,
                src_tif=current_f,
                out_prefix=tile_prefix,
            )
            
        # 3) add offsets
        if add_offsets:
            # extract local bounds from the ort_f
            #ullr_f = current_f.replace(".tif", f".bounds.offsets_x{x_offset}_y{y_offset}.txt")
            ullr_f =f"{tile_prefix}.globalcoord.ullr.txt"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Compute bounds + apply offsets")
            cmds.append(f"cartloader image_extract_ullr --input {current_f} --output {ullr_f} --x-offset {x_offset} --y-offset {y_offset}")
            mm.add_target(ullr_f, [current_f], cmds)
            # read the bounds
            gcoord_f = f"{tile_prefix}.globalcoord.tif"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Georeference (offset ULLR)")
            # append a bash command to get the ullr from the first row in ullr_f
            cmds.append(f"ULLR=$(head -n 1 {ullr_f}) && \\")
            cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr $ULLR {current_f} {gcoord_f}")
            mm.add_target(gcoord_f, [current_f, ullr_f], cmds)
            current_f = gcoord_f
        
        # 4) crop
        if args.crop_tile_by_minmax:
            # read the minmax
            xmin, ymin, xmax, ymax = idx2minmax.get((row, col), (0, 0, 0, 0))
            assert not (xmin == 0 and ymin == 0 and xmax == 0 and ymax == 0), f"Tile (row {row}, col {col}): min/max bounds not found (--in-minmax)"
            
            ullr_f = f"{tile_prefix}.crop.ullr.txt"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Compute crop bounds")
            cmds.append(f"cartloader image_extract_ullr --input {current_f} --output {ullr_f} --crop-by-minmax --minmax {xmin},{xmax},{ymin},{ymax}")
            mm.add_target(ullr_f, [current_f], cmds)
            
            # crop the image
            crop_f = f"{tile_prefix}.crop.tif"
            cmds = cmd_separator([], f"Tile (row {row}, col {col}): Crop by min/max")
            cmds.append(f"ULLR=$(head -n 1 {ullr_f}) && \\")
            cmds.append(f"{args.gdal_translate} -of GTiff -projwin $ULLR {current_f} {crop_f}")
            mm.add_target(crop_f, [current_f, ullr_f], cmds)
            current_f = crop_f

        tiles.append(current_f)

    # 3) stitch
    # gdalbuildvrt + gdal_translate
    cmds= cmd_separator([], f"Stitch tiles")
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
        cmds = cmd_separator([], f"Downsize image")
        downsize_percentage = f"{args.downsize_prop * 100}%"
        cmds.append(f"gdal_translate -outsize {downsize_percentage} {downsize_percentage} -r cubic {args.output} {downsize_f}")
        mm.add_target(downsize_f, [args.output], cmds)

    # write makefile
    if len(mm.targets) == 0:
        logging.error("There is no target to run. Please make sure that at least one run option was turned on")
        sys.exit(1)
    
    make_f = os.path.join(out_dir, args.makefn)
    mm.write_makefile(make_f)
    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)

if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
