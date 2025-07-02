import sys, os, gzip, argparse, subprocess
import pandas as pd
import numpy as np
import tifffile

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, add_param_to_cmd
from cartloader.utils.image_helper import orient2axisorder, update_orient

# get the current path
current_path = os.path.realpath(__file__)
cartloader_dir=os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
gdal_get_size_script = os.path.join(cartloader_dir, 'cartloader', "utils", "gdal_get_size.sh")

aux_image_arg={
    "transform": ["page", "level", "series", "upper_thres_quantile", "upper_thres_intensity", "lower_thres_quantile", "lower_thres_intensity", "colorize"]
    }

def check_ome_tiff(tiff_path):
    assert os.path.exists(tiff_path), f"The file '{tiff_path}' does not exist."
    with tifffile.TiffFile(tiff_path) as tif:
        if tif.is_ome:
            return True
        else:
            return False

def get_mono(args):
    # get the mono information
    with tifffile.TiffFile(args.in_fig) as tif:
        n_pages = len(tif.pages)
        if args.page is None:
            if n_pages > 1:
                raise ValueError("In --transform, multiple pages detected. Please specify the page number to extract the image from")
            elif n_pages == 1:
                args.page = 0
            else:
                raise ValueError("In --transform, no pages detected in the OME-TIFF file")
        page = tif.series[args.series].levels[args.level].pages[args.page]
        if len(page.shape) == 3:
            if page.shape[2] != 3:
                raise ValueError("In --transform, the colored image is not in RGB format")
            args.mono = False
        else:
            args.mono = True
    if args.colorize is not None:
        assert args.mono is True, "In --transform, the colorize option is only available for black-and-white images"
        args.mono = False
    return args.mono

def cmds_for_dimensions(geotif_f, dim_f):
    cmds = cmd_separator([], f"Extract dimensions from: {geotif_f}")
    dim_f = os.path.splitext(geotif_f)[0] + ".dim.tsv"
    print(geotif_f)
    print(dim_f)
    cmds.append(f"{gdal_get_size_script} {geotif_f} {dim_f}")
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
    parser = argparse.ArgumentParser(prog=f"cartloader run_fig2pmtiles", description="Convert a figure to pmtiles")

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it (default: False)')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default=None, help='The file name of the Makefile to generate (default: {out-prefix}.mk)')

    cmd_params = parser.add_argument_group("Commands", "Commands to run. If all functions were enabled, the steps will be executed in this order: transform, georeference, rotate, flip-vertical/horizontal, geotif2mbtiles, mbtiles2pmtiles, update-catalog")
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run main commands (geotif2mbtiles, mbtiles2pmtiles, upload-aws, update-yaml)')
    cmd_params.add_argument('--geotif2mbtiles', action='store_true', default=False, help='Convert a geotiff file to mbtiles')
    cmd_params.add_argument('--mbtiles2pmtiles', action='store_true', default=False, help='Convert mbtiles to pmtiles')
    #cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload the new pmtiles to the AWS S3') # Stop supporting this option
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Update the catalog.yaml file with the new pmtiles and upload to AWS')
    cmd_params.add_argument('--transform', action='store_true', default=False, help='Plus function. Convert Aligned Histology from 10x Xenium or Vizgen MERSCOPE. If enabled, the user must provide those transform arguments')
    cmd_params.add_argument('--georeference', action='store_true', default=False, help='Plus function. Create a geotiff file from PNG or TIF file. If enabled, the user must provide georeferenced bounds using --in-tsv or --in-bounds.')
    cmd_params.add_argument('--rotate', type=str, default=None, choices=["90", "180", "270"],  help='Plus function. Rotate the image by 90, 180, or 270 degrees clockwise. Rotate precedes flip.')
    cmd_params.add_argument('--flip-vertical', action='store_true', default=False, help='Plus function. Flip the image vertically (flipped along Y-axis). Rotate precedes flip.')
    cmd_params.add_argument('--flip-horizontal', action='store_true', default=False, help='Plus function. Flip the image horizontally (flipped along X-axis). Rotate precedes flip.')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Define the input file according to the user's needs.")
    inout_params.add_argument('--in-fig', type=str, help='The input figure file (PNG or TIF) to be converted to pmTiles')
    inout_params.add_argument('--out-prefix', required=True, type=str, help='The output prefix. A new directory will be created if needed')
    inout_params.add_argument('--in-tsv', type=str, default=None, help='If --georeference is required without --transform, use the *.pixel.sorted.tsv.gz from run_ficture to provide georeferenced bounds.')
    inout_params.add_argument('--in-bounds', type=str, default=None, help='If --georeference is required without --transform, provide the bounds in the format of "<ulx>,<uly>,<lrx>,<lry>", which represents upper-left X, upper-left Y, lower-right X, lower-right Y.')

    key_params = parser.add_argument_group("Key parameters")
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help='For --georeference and --geotif2mbtiles, define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--mono', action='store_true', default=False, help='For --rotate, --flip-vertical/horizontal, and --geotif2mbtiles without --transform, define if the input image is black-and-white and single-banded. (default: False)')
    key_params.add_argument('--rgba', action='store_true', default=False, help='RGBA, 4-banded image (default: False)')
    key_params.add_argument('--resample', type=str, default='cubic', help='For --geotif2mbtiles and --mbtiles2pmtiles, define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    key_params.add_argument('--blocksize', type=int, default='512', help='For --geotif2mbtiles and --mbtiles2pmtiles, define the blocksize (default: 512)')

    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')
    env_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate files')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters for steps, and tools")
    aux_params.add_argument('--catalog-yaml', type=str, default=None, help='For --update-catalog, define the catalog yaml file to update (default: catalog.yaml in the output directory specified by --out-prefix)')
    aux_params.add_argument('--basemap-key', type=str, default=None, help='For --update-catalog, The information to use for updating the basemap in the catalog.yaml file. An example for a histology file: <type>:<hist_id>.')
    aux_params.add_argument('--transform-csv', type=str, help='For --transform, (Vizgen only) CSV file containing conversion parameters (typically micron_to_mosaic_pixel_transform.csv)')
    aux_params.add_argument("--page", type=int, help='For --transform, for 3D (X/Y/Z) OME file, specify the z-value to extract the image from')
    aux_params.add_argument("--level", type=int, default=0, help='For --transform, specify the level to extract from the OME-TIFF file')
    aux_params.add_argument("--series", type=int, default=0, help='For --transform, specify the index of series to extract from the OME-TIFF file')
    aux_params.add_argument("--upper-thres-quantile", type=float, default=None, help='For --transform, quantile-based capped value for rescaling the image. Cannot be used with --upper-thres-intensity')
    aux_params.add_argument("--upper-thres-intensity", type=float, default=None, help='For --transform, intensity-based capped value for rescaling the image. Cannot be used with --upper-thres-quantile')
    aux_params.add_argument("--lower-thres-quantile", type=float, default=None, help='For --transform, quantile-based floored value for rescaling the image. Cannot be used with --lower-thres-intensity')
    aux_params.add_argument("--lower-thres-intensity", type=float, default=None, help='For --transform, intensity-based floored value for rescaling the image. Cannot be used with --lower-thres-quantile')
    aux_params.add_argument("--colorize", type=str, help='For --transform, colorize the black-and-white image using a specific RGB code as a max value (does not work with RGB images)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_fig2pmtiles(_args):

    args=parse_arguments(_args)

    if args.main:
        args.geotif2mbtiles = True
        args.mbtiles2pmtiles = True
        #args.upload_aws = True
        args.update_catalog = True

    if args.update_catalog and args.basemap_key is None:
        raise ValueError("Please provide the basemap key using --basemap-key for updating the catalog.yaml file")

    scheck_app(args.pmtiles)
    scheck_app(args.gdal_translate)
    scheck_app(args.gdaladdo)

    # files
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir != "":
        os.makedirs(out_dir, exist_ok=True)
    # - use in_fig as the input figure for all steps, which will be constantly updated
    assert args.in_fig is not None and os.path.exists(args.in_fig), "Please provide a valid input figure file using --in-figure"
    # - mbtiles
    mbtile_f=f"{args.out_prefix}.pmtiles.mbtiles"
    mbtile_flag = f"{mbtile_f}.done"
    mbtile_f_ann=f"{args.out_prefix}.pmtiles.{args.resample}.mbtiles"
    # - pmtiles
    pmtiles_f=f"{args.out_prefix}.pmtiles"
    catalog_f = f"{out_dir}/catalog.yaml" if args.catalog_yaml is None else args.catalog_yaml

    # start mm
    mm = minimake()
    
    # 0. Check the orientation, and update it to the equivalent transformations
    if args.flip_vertical or args.flip_horizontal or args.rotate is not None:
        args.rotate, args.flip_vertical, args.flip_horizontal = update_orient (args.rotate, args.flip_vertical, args.flip_horizontal, args.in_fig)
        
    # 1. Transform the figure
    if args.transform:
        # update the mono information from transform
        args.mono = get_mono(args)
        #assert check_ome_tiff(args.in_fig), "When --transform is enabled, the input figure must be an OME-TIFF file"
        transform_prefix=f"{args.out_prefix}.transform"
        transform_f=f"{transform_prefix}.png"
        cmds = cmd_separator([], f"Transforming {args.in_fig} to {transform_f}")
        cmd= " ".join([
            "cartloader transform_aligned_histology",
            f"--tif {args.in_fig}", 
            f"--out-prefix {transform_prefix}",
            f"--csv {args.transform_csv}" if args.transform_csv is not None else "",
            f"--rotate-clockwise" if args.rotate == "90" else "",
            f"--rotate-counter" if args.rotate == "270" else "",
            f"--flip-vertical" if args.flip_vertical else "",
            f"--flip-horizontal" if args.flip_horizontal else "",
            "--skip-pmtiles",
            "--write-bounds"
        ])
        cmd = add_param_to_cmd(cmd, args, aux_image_arg["transform"])
        cmds.append(cmd)
        # >1 output: transform_f, transform_prefix.bounds 
        cmds.append(f"[ -f {transform_f} ] && [ -f {transform_prefix}.bounds.csv ] && touch {transform_prefix}.done")
        mm.add_target(f"{transform_prefix}.done", [args.in_fig], cmds)
    else:
        transform_f = args.in_fig

    # 2. Create a geotiff file from PNG or TIF file
    if args.georeference:
        georef_f =f"{args.out_prefix}.georef.tif"
        cmds = cmd_separator([], f"Geo-referencing {transform_f} to {georef_f}")

        if args.transform:
            # add a bash cmd to read ulx, uly, lrx, lry from the bounds file, which is one line file with the format of "ulx,uly,lrx,lry"
            cmds.append(f"bounds=$(sed 's/,/ /g' {transform_prefix}.bounds.csv) ")
            cmd = " ".join([
                "gdal_translate",
                "-of GTiff",
                "-a_srs", args.srs,
                "-a_ullr", "$bounds",
                transform_f, 
                georef_f
            ])
            #cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr $bounds {transform_f} {georef_f}")
            mm.add_target(georef_f, [f"{transform_prefix}.done"], cmds)
            # update the flip/rotate information since the transform will be carried out in georeferencing.
            args.rotate=None
            args.flip_vertical=False
            args.flip_horizontal=False
        else:
            if  args.in_bounds is not None:
                ulx,uly,lrx,lry = args.in_bounds.split(",")
            elif args.in_tsv is not None:
                with gzip.open(args.in_tsv, 'rt') as f:
                    for i in range(3):
                        line = f.readline()
                ann2val = {x.split("=")[0]:x.split("=")[1] for x in line.strip().replace("##", "").split(";")}
                ulx = float(ann2val["OFFSET_X"]) # Upper left X-coordinate
                uly = float(ann2val["OFFSET_Y"]) # Upper left Y-coordinate
                lrx = float(ann2val["SIZE_X"])+1+ulx # Lower Right X-coordinate
                lry = float(ann2val["SIZE_Y"])+1+uly # Lower Right Y-coordinate 
            else:
                raise ValueError("Please provide either --in-bounds or --in-tsv to georeference the figure")

            cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr {ulx} {uly} {lrx} {lry} {transform_f} {georef_f}")
            mm.add_target(georef_f, [transform_f], cmds)
    else:
        georef_f = transform_f

    # 3. Orientation
    if args.flip_vertical or args.flip_horizontal or args.rotate is not None:
        # 0) Extract dimensions
        dim_f=georef_f.replace(".tif","") + ".dim.tsv"
        cmds=cmds_for_dimensions(georef_f, dim_f)
        mm.add_target(dim_f, [georef_f], cmds)
    
        # 1) Output FILENAME
        ort_suffix=get_orientation_suffix(args.rotate, args.flip_vertical, args.flip_horizontal)
        ort_f = f"{args.out_prefix}.{ort_suffix}.tif"

        # 2) Orientation
        cmds =  cmds_for_orientation(georef_f, dim_f, ort_f, args.rotate, args.flip_vertical, args.flip_horizontal, args.mono, args.rgba)
        mm.add_target(ort_f, [georef_f, dim_f], cmds)
    else:
        ort_f = georef_f

    # 4. Convert a raster image (GeoTIFF) into map tiles (MBTiles format)
    # * MBTiles expects top-left origin. If the input geotif use lower left was origin, GDAL will flip the origin to top left. 
    # * The origin could be find by running `gdalinfo <filename>` and check the `Pixel Size=(<x>,<y>)` field. If Pixel Y is negative, image is stored top-down. 
    # * Use a flag file instead of mbtiles_f to indicate the completion of the mbtiles conversion - to avoid silently failed or interrupted conversions.
    if args.geotif2mbtiles:
        cmds = cmd_separator([], f"Converting from geotif to mbtiles: {ort_f}")
        partial_db = mbtile_f.replace('.mbtiles', '.partial_tiles.db')
        journal_db = f"{mbtile_f}-journal"

        create_mbtile_flag(mbtile_flag, mbtile_f, partial_db, journal_db)

        ## Add cleanup step to prevent the case that the previous run was interrupted and left behind a .mbtiles file with a journal file or a partial_tiles.db file. Such unfinished jobs will not be availble to detect by the makefile.
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
            "-scale",                                       # Automatically scale pixel values to 0-255
            "-of", "mbtiles",
            "-a_srs", args.srs,                             # default to EPSG:3857. Assign the target projection: Web Mercator (used by web maps)
            ort_f, 
            mbtile_f
        ])
        cmds.append(cmd)
        
        validation_cmd = f" [ -f {mbtile_f} ]  && [ ! -f {journal_db} ] && [ ! -f {partial_db} ] && touch {mbtile_flag}"
        cmds.append(validation_cmd)
        mm.add_target(mbtile_flag, [ort_f], cmds)
    
    # 5. Convert mbtiles to pmtiles
    if args.mbtiles2pmtiles:
        cmds = cmd_separator([], f"Resampling mbtiles and converting to pmtiles: {mbtile_f}")
        cmds.append(f"cp {mbtile_f} {mbtile_f_ann}")
        cmds.append(f"'{args.gdaladdo}' {mbtile_f_ann} -r {args.resample} 2 4 8 16 32 64 128 256") # Build internal overviews. 
        cmds.append(f"'{args.pmtiles}' convert --force {mbtile_f_ann} {pmtiles_f}")
        mm.add_target(pmtiles_f, [mbtile_flag], cmds)

    # 6. Update the catalog.yaml file with the new pmtiles and upload to AWS
    if args.update_catalog:
        cmds = cmd_separator([], f"Updating yaml for pmtiles: {pmtiles_f}")
        pmtiles_name=os.path.basename(pmtiles_f)
        pmtiles_dir=os.path.dirname(pmtiles_f)
        cmd= " ".join([
            "cartloader write_catalog_for_assets",
            f"--out-catalog {catalog_f}",
            f"--basemap {args.basemap_key}:{pmtiles_name}",
            f"--basemap-dir {pmtiles_dir}"
        ])
        cmds.append(cmd)
        #cmds.append(f"cartloader update_catalog_for_basemap --in-yaml {catalog_f} --basemap {args.basemap_key}:{pmtiles_name} --basemap-dir {pmtiles_dir} --overwrite")
        cmds.append(f"touch {pmtiles_f}.yaml.done")
        mm.add_target(f"{pmtiles_f}.yaml.done", [pmtiles_f, catalog_f], cmds)

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
    make_f=os.path.join(out_dir, args.makefn) if args.makefn is not None else f"{args.out_prefix}.mk"
    mm.write_makefile(make_f)

    if args.dry_run:
        dry_cmd=f"make -f {make_f} -n {'-B' if args.restart else ''} "
        os.system(dry_cmd)
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs}")
    else:
        exe_cmd=f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}"
        result = subprocess.run(exe_cmd, shell=True)
        if result.returncode != 0:
            print(f"Error in executing: {exe_cmd}")
            sys.exit(1)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])