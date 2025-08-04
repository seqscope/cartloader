import sys, os, gzip, argparse, subprocess, inspect
import pandas as pd
import numpy as np
import tifffile

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, add_param_to_cmd, write_dict_to_file, load_file_to_dict, scheck_actions, execute_makefile
from cartloader.utils.image_helper import orient2axisorder, update_orient
from cartloader.scripts.image_png2pmtiles import get_orientation_suffix

# get the current path
current_path = os.path.realpath(__file__)
cartloader_dir=os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
gdal_get_size_script = os.path.join(cartloader_dir, 'cartloader', "utils", "gdal_get_size.sh")

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="Convert a figure (OME-TIFF, TIFF, or PNG) to pmtiles")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default=None, help='The file name of the Makefile to generate (default: {out-prefix}.mk)')

    cmd_params = parser.add_argument_group("Commands", "Commands to run. In the order of: ome2png, georeference, rotate/flip, png2pmtiles")
    cmd_params.add_argument('--ome2png', action='store_true', default=False, help='Convert OME-TIFF to PNG (Typically used for Vizgen and Xenium). If enabled, this will automatically set --georeference with bounds information from OME-TIFF.')
    cmd_params.add_argument('--png2pmtiles', action='store_true', default=False, help='Convert PNG to pmtiles')
    cmd_params.add_argument('--georeference', action='store_true', default=False, help='Plus function. Create a geotiff file from PNG or TIF file. If enabled, the user must provide georeferenced bounds using --in-tsv or --in-bounds.')
    cmd_params.add_argument('--rotate', type=str, default=None, choices=["90", "180", "270"],  help='Plus function. Rotate the image by 90, 180, or 270 degrees clockwise. Rotate precedes flip.')
    cmd_params.add_argument('--flip-vertical', action='store_true', default=False, help='Plus function. Flip the image vertically (flipped along Y-axis). Rotate precedes flip.')
    cmd_params.add_argument('--flip-horizontal', action='store_true', default=False, help='Plus function. Flip the image horizontally (flipped along X-axis). Rotate precedes flip.')
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Plus function. Flip the image horizontally (flipped along X-axis). Rotate precedes flip.')

    inout_params = parser.add_argument_group("Input/Output Parameters", """
                                             Two ways to define the input and output: 
                                             1) use the --in-json and --fig-id to locate input image in the JSON file.
                                             2) use --in-img to provide the path to input image.
                                             """)
    inout_params.add_argument('--out-dir', type=str, default=None, required=True, help='The output directory.')
    inout_params.add_argument('--in-json', type=str, default=None, help='Input JSON. It should map figure IDs to image paths. Each key is a figure ID, and each value is the corresponding file path. If this is provided, --in-img will be ignored.')
    inout_params.add_argument('--in-img', type=str, help='The input image file (PNG or TIF) to be converted to pmTiles')
    inout_params.add_argument('--img-id', type=str, default=None, required=True, help='Image ID. This will be used as the output file name prefix. Also, if --in-json is provided, this will be used to locate the image in the JSON file.')
    inout_params.add_argument('--catalog-yaml', type=str, default=None, help='For --update-catalog, define the catalog yaml file to update (default: <out_dir>/catalog.yaml)')

    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    env_params.add_argument('--gdalinfo', type=str, default=f"gdalinfo", help='Path to gdalinfo binary')

    aux_params1 = parser.add_argument_group("Auxiliary parameters for --ome2png")
    aux_params1.add_argument('--micron2pixel-csv', type=str, help='(Vizgen only) CSV file containing transformation parameters from microns to mosaic pixels. (typically micron_to_mosaic_pixel_transform.csv)')
    aux_params1.add_argument("--page", type=int, help='For 3D (X/Y/Z) OME file, specify the z-value to extract the image from')
    aux_params1.add_argument("--level", type=int, help='Level to extract from the OME-TIFF file')
    aux_params1.add_argument("--series", type=int, help='Index of series to extract from the OME-TIFF file')
    aux_params1.add_argument("--upper-thres-quantile", type=float, default=None, help='Quantile-based capped value for rescaling the image. Cannot be used with --upper-thres-intensity')
    aux_params1.add_argument("--upper-thres-intensity", type=float, default=None, help='Intensity-based capped value for rescaling the image. Cannot be used with --upper-thres-quantile')
    aux_params1.add_argument("--lower-thres-quantile", type=float, default=None, help='Quantile-based floored value for rescaling the image. Cannot be used with --lower-thres-intensity')
    aux_params1.add_argument("--lower-thres-intensity", type=float, default=None, help='Intensity-based floored value for rescaling the image. Cannot be used with --lower-thres-quantile')
    aux_params1.add_argument("--transparent-below", type=int, default=None, help='Set pixels below this value to transparent (0-255)')
    aux_params1.add_argument("--colorize", type=str, help='Colorize the black-and-white image using a specific RGB code as a max value (does not work with RGB images)')
    aux_params1.add_argument('--high-memory', action='store_true', default=False)

    aux_params2 = parser.add_argument_group("Auxiliary parameters for for --png2pmtiles")
    aux_params2.add_argument('--georef-pixel-tsv', type=str, default=None, help='If --georeference is required without --ome2png, use the *.pixel.sorted.tsv.gz from run_ficture to provide georeferenced bounds.')
    aux_params2.add_argument('--georef-bounds-tsv', type=str, default=None, help='If --georeference is required without --ome2png, use a tsv file with one line of <ulx>,<uly>,<lrx>,<lry> to provide georeferenced bounds.')
    aux_params2.add_argument('--georef-bounds', type=str, default=None, help='If --georeference is required without --ome2png, provide the bounds in the format of "<ulx>,<uly>,<lrx>,<lry>", which represents upper-left X, upper-left Y, lower-right X, lower-right Y.')
    aux_params2.add_argument('--srs', type=str, default='EPSG:3857', help='For --png2pmtiles, define the spatial reference system (default: EPSG:3857)')
    aux_params2.add_argument('--mono', action='store_true', default=False, help='For --png2pmtiles without --ome2png, define if the input image is black-and-white and single-banded manually. (default: False)')
    aux_params2.add_argument('--rgba', action='store_true', default=False, help='For --png2pmtiles without --ome2png, define if the image is RGBA, i.e., 4-banded, image manually (default: False)')
    aux_params2.add_argument('--resample', type=str, default='cubic', help='For --png2pmtiles, define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    aux_params2.add_argument('--blocksize', type=int, default='512', help='For --png2pmtiles, define the blocksize (default: 512)')


    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

aux_image_arg={
    "ome2png": ["page", "level", "series", "upper_thres_quantile", "upper_thres_intensity", "lower_thres_quantile", "lower_thres_intensity", "transparent_below", "colorize", "high_memory"],
    "png2pmtiles": ["srs", "mono", "rgba", "resample", "blocksize", "pmtiles", "gdaladdo"],
    "georeference": ["georef_pixel_tsv", "georef_bounds_tsv", "georef_bounds", "srs"],
    "orientate": ["gdalinfo"]
}


def get_mono(args):
    # get the mono information
    with tifffile.TiffFile(args.in_img) as tif:
        n_pages = len(tif.pages)
        if args.page is None:
            if n_pages > 1:
                raise ValueError("In --ome2png, multiple pages detected. Please specify the page number to extract the image from")
            elif n_pages == 1:
                args.page = 0
            else:
                raise ValueError("In --ome2png, no pages detected in the OME-TIFF file")
        page = tif.series[args.series].levels[args.level].pages[args.page]
        if len(page.shape) == 3:
            if page.shape[2] != 3:
                raise ValueError("In --ome2png, the colored image is not in RGB format")
            args.mono = False
        else:
            args.mono = True
    if args.colorize is not None:
        assert args.mono is True, "In --ome2png, the colorize option is only available for black-and-white images"
        args.mono = False
    return args.mono


def import_image(_args):
    args=parse_arguments(_args)

    # actions
    scheck_actions(args, ["--ome2png", "--png2pmtiles", "--georeference", "--rotate", "--flip-vertical", "--flip-horizontal"], context="actions")

    # output dir
    os.makedirs(args.out_dir, exist_ok=True)
    img_prefix =  os.path.join(args.out_dir, args.img_id)

    # input files
    # read in_json if provided
    assert args.img_id, f"Provide an ID for the image using --img-id"
    if args.in_json is not None:
        assert os.path.exists(args.in_json), f"The input json file doesn't exist: {args.in_json}"
        xenium_ranger_data=load_file_to_dict(args.in_json)
        # only keep the key:value pairs that have a key in ["DAPI", "BOUNDARY_IMG", "INTERIOR_RNA_IMG", "INTERIOR_PROTEIN_IMG", "DAPI_3D", "DAPI_MIP"]
        args.in_img = xenium_ranger_data[args.img_id]

    assert args.in_img is not None and os.path.exists(args.in_img), "Please provide a valid input figure file using --in-img"
    
    # start mm
    mm = minimake()
    
    # 0. Check the orientation, and update it to the equivalent transformations
    if args.flip_vertical or args.flip_horizontal or args.rotate is not None:
        args.rotate, args.flip_vertical, args.flip_horizontal = update_orient (args.rotate, args.flip_vertical, args.flip_horizontal)
        
    # 1. Transform the figure
    if args.ome2png:
        transform_prefix = f"{img_prefix}.transform"
        transform_bounds_tsv = f"{transform_prefix}.bounds.csv"   
        transform_f = f"{transform_prefix}.png"
        color_mode=f"{transform_prefix}.color.csv"

        cmds = cmd_separator([], f"Converting OME TIFF ({args.in_img}) to PNG ({transform_f})")
        cmd= " ".join([
            "cartloader image_ome2png",
            f"--tif {args.in_img}", 
            f"--csv {args.micron2pixel_csv}" if args.micron2pixel_csv is not None else "",
            f"--out-prefix {transform_prefix}",
            f"--rotate-clockwise" if args.rotate == "90" else "", # Note image_ome2png has different options than the 90,180,270
            f"--rotate-counter" if args.rotate == "270" else "",
            f"--flip-vertical" if args.flip_vertical else "",
            f"--flip-horizontal" if args.flip_horizontal else "",
            f"--write-color-mode"
        ])
        cmd = add_param_to_cmd(cmd, args, aux_image_arg["ome2png"])
        cmds.append(cmd)
        # >1 output: transform_f, transform_prefix.bounds 
        cmds.append(f"[ -f {transform_f} ] && [ -f {transform_bounds_tsv} ] && touch {transform_prefix}.done")
        mm.add_target(f"{transform_prefix}.done", [args.in_img], cmds)
        
        # update for georeference and bounds
        args.georeference = True
        assert not args.georef_bounds and not args.georef_pixel_tsv and not args.georef_bounds_tsv, f"Since --ome2png, skip --georef-pixel-tsv, --georef-bounds, and --georef-bounds-tsv. The georeferenced bounds will be automatically extract from the OME TIFF file"
        args.georef_bounds_tsv =  transform_bounds_tsv

        # update the flip/rotate information and georeferencing
        args.rotate=None
        args.flip_vertical=False
        args.flip_horizontal=False
        
        # update input for the next step
        flat_img = transform_f
    else:
        flat_img = args.in_img
        color_mode=None

    # 2. Use image png2pmtiles
    pmtiles_f = os.path.join(args.out_dir, f"{args.img_id}.pmtiles")
    if args.png2pmtiles:
        if args.ome2png:
            prereq = [f"{transform_prefix}.done"]
        else:
            prereq = [flat_img]
        cmds = cmd_separator([], f"Converting PNG ({flat_img}) to PMTiles ({pmtiles_f})")
        cmd = " ".join([
            "cartloader image_png2pmtiles",
            # actions
            "--georeference" if args.georeference else "",
            f"--rotate {args.rotate}" if args.rotate else "",
            "--flip-vertical" if args.flip_vertical else "",
            "--flip-horizontal" if args.flip_horizontal else "", 
            "--geotif2mbtiles",
            "--mbtiles2pmtiles", # this step will write down a asset JSON file
            # in
            f"--in-img {flat_img}",
            # out
            f"--out-prefix {img_prefix}",
            # params
            f"--color-mode-record {color_mode}" if color_mode else "",
            f"--mono {args.mono}" if args.mono and not color_mode else "",
            f"--rgba {args.rgba}" if args.rgba and not color_mode else "",
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
        ])
        cmd = add_param_to_cmd(cmd, args, list(set(aux_image_arg["png2pmtiles"] + aux_image_arg["georeference"])))
        cmd = add_param_to_cmd(cmd, args, ["restart", "n_jobs"])
        cmds.append(cmd)
        mm.add_target(pmtiles_f, prereq, cmds)
    elif args.georeference or args.flip_vertical or args.flip_horizontal or args.rotate is not None:
        # define actions
        png2pmtiles_actions=[]
        if args.georeference:
            png2pmtiles_actions.append("georeferencing")
        if args.flip_vertical or args.flip_horizontal or args.rotate is not None:
            png2pmtiles_actions.append("orientation")

        # define the output from the last step
        if "georeferencing" in png2pmtiles_actions and "orientation" not in png2pmtiles_actions:
            cond_out =f"{img_prefix}.georef.tif"
        else:
            ort_suffix=get_orientation_suffix(args.rotate, args.flip_vertical, args.flip_horizontal)
            cond_out = f"{img_prefix}.{ort_suffix}.tif"
        
        cmds = cmd_separator([], f"Processing PNG ({flat_img}): {", ".join(png2pmtiles_actions)}")
        cmd = " ".join([
            "cartloader image_png2pmtiles",
            "--georeference" if args.georeference else "",
            f"--rotate {args.rotate}" if args.rotate else "",
            "--flip-vertical" if args.flip_vertical else "",
            "--flip-horizontal" if args.flip_horizontal else "", 
            f"--in-img {flat_img}",
            f"--out-prefix {img_prefix}",
            f"--color-mode-record {color_mode}" if color_mode else "",
            f"--mono {args.mono}" if args.mono and not color_mode else "",
            f"--rgba {args.rgba}" if args.rgba and not color_mode else "",
            f"--gdal_translate {args.gdal_translate}" if args.gdal_translate else "",
        ])
        cmd = add_param_to_cmd(cmd, args, list(set(aux_image_arg["orientate"] + aux_image_arg["georeference"])))
        cmds.append(cmd)
        mm.add_target(cond_out, prereq, cmds)
    
    # 3. Update the catalog.yaml file with the new pmtiles and upload to AWS
    if args.update_catalog:
        if args.catalog_yaml is None:
            args.catalog_yaml = os.path.join(args.out_dir, "catalog.yaml")
        cmds = cmd_separator([], f"Updating yaml for pmtiles: {pmtiles_f}")        
        cmd= " ".join([
            "cartloader write_catalog_for_assets",
            f"--out-catalog {args.catalog_yaml}",
            f"--basemap {args.img_id}:{args.img_id}.pmtiles",
            f"--basemap-dir {args.out_dir}"
        ])
        cmds.append(cmd)
        cmds.append(f"touch {pmtiles_f}.yaml.done")
        mm.add_target(f"{pmtiles_f}.yaml.done", [pmtiles_f, args.catalog_yaml], cmds)

    # 7. Upload new PMtiles to AWS
    # if args.upload_aws:
    #     assert args.aws_dir is not None, "Please provide the AWS S3 bucket path using --aws-bucket"
    #     cmds = cmd_separator([], f"Uploading pmtiles to AWS: {geotif_f}")
    #     pmtiles_fn=os.path.basename(pmtiles_f)
    #     cmds.append(f"aws s3 cp {pmtiles_f} {args.aws_dir}/{pmtiles_fn}")
    #     cmds.append(f"aws s3 cp {catalog_f} {args.aws_dir}/{catalog_f}")
    #     cmds.append(f"touch {pmtiles_f}.aws.done")
    #     mm.add_target(f"{pmtiles_f}.aws.done", [pmtiles_f], cmds)

    ## write makefile
    make_f=os.path.join(args.out_dir, args.makefn) if args.makefn is not None else f"{img_prefix}.mk"
    mm.write_makefile(make_f)

    execute_makefile(make_f, dry_run=args.dry_run, restart=args.restart, n_jobs=args.n_jobs)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])