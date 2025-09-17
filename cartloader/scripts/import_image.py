import sys, os, gzip, argparse, subprocess, inspect
import pandas as pd
import numpy as np
import tifffile

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, add_param_to_cmd, write_dict_to_file, load_file_to_dict, scheck_actions, execute_makefile
from cartloader.utils.image_helper import orient2axisorder, update_orient
from cartloader.scripts.image_png2pmtiles import get_orientation_suffix


def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Import an image (OME-TIFF/TIFF/PNG) and produce PMTiles; optional georeferencing/orientation and catalog update"
    )

    run_params = parser.add_argument_group("Run Options", "Execution controls for generating and running the Makefile")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Generate the Makefile but do not execute it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Ignore existing outputs and re-run all steps')
    run_params.add_argument('--makefn', type=str, default=None, help='Name of the generated Makefile (default: {out-prefix}.mk)')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of parallel jobs to run (default: 1)')

    cmd_params = parser.add_argument_group("Commands", "Actions to apply (order: ome2png → georeference/orientation → png2pmtiles)")
    cmd_params.add_argument('--ome2png', action='store_true', default=False, help='Convert OME-TIFF to PNG (e.g., Vizgen/Xenium). Also enables --georeference using bounds read from the OME-TIFF')
    cmd_params.add_argument('--png2pmtiles', action='store_true', default=False, help='Generate PMTiles from a PNG (and intermediate GeoTIFF)')
    cmd_params.add_argument('--georeference', action='store_true', default=False, help='Georeference the image via bounds; required unless bounds come from --ome2png; provide bounds with --georef-*')
    cmd_params.add_argument('--rotate', type=str, default=None, choices=["90", "180", "270"],  help='Rotate clockwise by 90/180/270 degrees (applied before flips)')
    cmd_params.add_argument('--flip-vertical', action='store_true', default=False, help='Flip vertically (around X axis); applied after rotation')
    cmd_params.add_argument('--flip-horizontal', action='store_true', default=False, help='Flip horizontally (around Y axis); applied after rotation')
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Update catalog.yaml with the generated PMTiles')

    inout_params = parser.add_argument_group(
        "Input/Output Parameters",
        'Two input modes: 1) JSON + --img-id (use --in-json to map image IDs to paths); 2) direct path via --in-img'
    )
    inout_params.add_argument('--out-dir', type=str, required=True, help='Output directory')
    inout_params.add_argument('--in-json', type=str, default=None, help='JSON/YAML mapping from image IDs to file paths; when set, --img-id selects the entry and --in-img is ignored')
    inout_params.add_argument('--in-img', type=str, help='Path to input image (PNG or OME-TIFF/TIFF)')
    inout_params.add_argument('--img-id', type=str, required=True, help='Image ID used as output filename prefix. When used with --in-json, selects the image to process')
    inout_params.add_argument('--catalog-yaml', type=str, default=None, help='Catalog YAML path (required if --update-catalog; default: <out-dir>/catalog.yaml)')

    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles (default: pmtiles)')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary (default: gdaladdo)')
    env_params.add_argument('--gdalinfo', type=str, default=f"gdalinfo", help='Path to gdalinfo binary (default: gdalinfo)')

    aux_params1 = parser.add_argument_group("Auxiliary parameters for --ome2png")
    aux_params1.add_argument('--micron2pixel-csv', type=str, help='CSV file containing transformation parameters from microns to mosaic pixels (platform: Vizgen; typical: micron_to_mosaic_pixel_transform.csv)')
    aux_params1.add_argument("--page", type=int, help='Z-slice index to extract from multi-page OME-TIFF (3D)')
    aux_params1.add_argument("--level", type=int, help='Resolution level index to extract from OME-TIFF')
    aux_params1.add_argument("--series", type=int, help='Series index to extract from OME-TIFF')
    aux_params1.add_argument("--upper-thres-quantile", type=float, default=None, help='Rescale cap by quantile; mutually exclusive with --upper-thres-intensity')
    aux_params1.add_argument("--upper-thres-intensity", type=float, default=None, help='Rescale cap by intensity; mutually exclusive with --upper-thres-quantile')
    aux_params1.add_argument("--lower-thres-quantile", type=float, default=None, help='Rescale floor by quantile; mutually exclusive with --lower-thres-intensity')
    aux_params1.add_argument("--lower-thres-intensity", type=float, default=None, help='Rescale floor by intensity; mutually exclusive with --lower-thres-quantile')
    aux_params1.add_argument("--transparent-below", type=int, default=None, help='Set pixels below this value to transparent (range: 0-255)')
    aux_params1.add_argument("--colorize", type=str, help='Colorize a mono image using an RGB hex or name. Cannot be used with RGB images.')
    aux_params1.add_argument('--high-memory', action='store_true', default=False, help='Use a memory-intensive path for better performance on large OME-TIFFs')

    aux_params2 = parser.add_argument_group("Auxiliary parameters for --png2pmtiles")
    aux_params2.add_argument('--srs', type=str, default='EPSG:3857', help='Spatial reference system identifier (default: EPSG:3857)')
    aux_params2.add_argument('--mono', action='store_true', default=False, help='Input is single-band (mono) PNG. Skip it if --ome2png is enabled.')
    aux_params2.add_argument('--rgba', action='store_true', default=False, help='Input is 4-band RGBA PNG. Skip it if --ome2png is enabled')
    aux_params2.add_argument('--resample', type=str, default='cubic', help='Resampling method: near, bilinear, cubic, etc (default: cubic)')
    aux_params2.add_argument('--blocksize', type=int, default='512', help='Block size in pixels for GDAL operations (default: 512)')

    aux_params3 = parser.add_argument_group("Auxiliary parameters for --georeference", "Pick one of the following three ways to provide georeferencing bounds")
    aux_params3.add_argument('--georef-pixel-tsv', type=str, default=None, help='Bounds source: *.pixel.sorted.tsv.gz from run_ficture2. Skip it if --ome2png is enabled')
    aux_params3.add_argument('--georef-bounds-tsv', type=str, default=None, help='Bounds source TSV with one line: <ulx>,<uly>,<lrx>,<lry>. Skip it if --ome2png is enabled')
    aux_params3.add_argument('--georef-bounds', type=str, default=None, help='Bounds string: "<ulx>,<uly>,<lrx>,<lry>". Skip it if --ome2png is enabled')

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


# def get_mono(args):
#     # get the mono information
#     with tifffile.TiffFile(args.in_img) as tif:
#         n_pages = len(tif.pages)
#         if args.page is None:
#             if n_pages > 1:
#                 raise ValueError("In --ome2png, multiple pages detected. Please specify the page number to extract the image from")
#             elif n_pages == 1:
#                 args.page = 0
#             else:
#                 raise ValueError("In --ome2png, no pages detected in the OME-TIFF file")
#         page = tif.series[args.series].levels[args.level].pages[args.page]
#         if len(page.shape) == 3:
#             if page.shape[2] != 3:
#                 raise ValueError("In --ome2png, the colored image is not in RGB format")
#             args.mono = False
#         else:
#             args.mono = True
#     if args.colorize is not None:
#         assert args.mono is True, "In --ome2png, the colorize option is only available for black-and-white images"
#         args.mono = False
#     return args.mono


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
        assert os.path.exists(args.in_json), f"File not found: {args.in_json} (--in-json)"
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
        
        cmds = cmd_separator([], f"Processing PNG ({flat_img}): {',' .join(png2pmtiles_actions)}")
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
