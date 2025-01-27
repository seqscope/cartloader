import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, re
import pandas as pd
import numpy as np

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger

# get the current path
current_path = os.path.realpath(__file__)
cartloader_dir=os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
gdal_get_size_script = os.path.join(cartloader_dir, 'cartloader', "utils", "gdal_get_size.sh")

def cmds_for_dimensions(geotif_f, dim_f):
    cmds = cmd_separator([], f"Extract dimensions from: {geotif_f}")
    dim_f = geotif_f.replace(".tif",".dim.tsv") 
    cmds.append(f"{gdal_get_size_script} {geotif_f} {dim_f}")
    return cmds

rotation_flip_map = {
    # rotate, flip_vertical, flip_horizontal: order_arg
    # no rotation,
    (None, False, True): "-1,2",     # Flip X-axis
    (None, True, False): "1,-2",     # Flip Y-axis
    (None, True, True): "-1,-2",     # Flip both axes

    # Rotate 90° clockwise
    ("90", False, False): "2,-1",    # Rotate 90° clockwise
    ("90", False, True): "-2,-1",    # Rotate 90°, then flip X
    ("90", True, False): "2,1",      # Rotate 90°, then flip Y
    ("90", True, True): "-2,1",      # Rotate 90°, then flip both

    # Rotate 180° clockwise
    ("180", False, False): "-1,-2",  # Rotate 180° clockwise
    ("180", False, True): "1,-2",    # Rotate 180°, then flip X
    ("180", True, False): "-1,2",    # Rotate 180°, then flip Y
    ("180", True, True): "1,2",      # Rotate 180°, then flip both (cancels out)

    # Rotate 270° clockwise
    ("270", False, False): "-2,1",   # Rotate 270° clockwise
    ("270", False, True): "2,1",     # Rotate 270°, then flip X
    ("270", True, False): "-2,-1",   # Rotate 270°, then flip Y
    ("270", True, True): "2,-1",     # Rotate 270°, then flip both
}

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    # repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader run_fig2pmtiles", description="Convert a figure to pmtiles")

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--main', action='store_true', default=False, help='Run main commands (geotif2mbtiles, mbtiles2pmtiles, upload-aws, update-yaml)')
    cmd_params.add_argument('--geotif2mbtiles', action='store_true', default=False, help='Convert a geotiff file to mbtiles')
    cmd_params.add_argument('--mbtiles2pmtiles', action='store_true', default=False, help='Convert mbtiles to pmtiles')
    #cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload the new pmtiles to the AWS S3') # Stop supporting this option
    cmd_params.add_argument('--update-catalog', action='store_true', default=False, help='Update the catalog.yaml file with the new pmtiles and upload to AWS')
    cmd_params.add_argument('--georeference', action='store_true', default=False, help='Plus function. Create a geotiff file from PNG or TIF file. If enabled, the user must provide georeferenced bounds using --in-tsv or --in-bounds.')
    cmd_params.add_argument('--rotate', type=str, default=None, choices=["90", "180", "270"],  help='Plus function. Rotate the image by 90, 180, or 270 degrees clockwise. If both rotate and flip options are provided, the rotation is applied first. ')
    cmd_params.add_argument('--flip-vertical', action='store_true', default=False, help='Plus function. Flip the image vertically (flipped along Y-axis). If both rotate and flip options are provided, the rotation is applied first. ')
    cmd_params.add_argument('--flip-horizontal', action='store_true', default=False, help='Plus function. Flip the image horizontally (flipped along X-axis). If both rotate and flip options are provided, the rotation is applied first.')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Define the input file according to the user's needs.")
    inout_params.add_argument('--in-fig', type=str, help='The input figure file (PNG or TIF) to be converted to pmTiles')
    inout_params.add_argument('--in-tsv', type=str, default=None, help='If --georeference is required, use the *.pixel.sorted.tsv.gz from run_ficture to provide georeferenced bounds.')
    inout_params.add_argument('--in-bounds', type=str, default=None, help='If --georeference is required, provide the bounds in the format of "<ulx>,<uly>,<lrx>,<lry>", which represents upper-left X, upper-left Y, lower-right X, lower-right Y.')
    inout_params.add_argument('--out-prefix', required=True, type=str, help='The output prefix. TNew directory will be created if needed')
    inout_params.add_argument('--basemap-key', type=str, default=None, help='The information to use for updating the basemap in the catalog.yaml file. An example for a histology file: <type>:<hist_id>.')

    key_params = parser.add_argument_group("Key parameters")
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help='For the georeference and geo2tiff steps, define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--mono', action='store_true', default=False, help='Black-and-white, single-banded image (default: False)')
    key_params.add_argument('--resample', type=str, default='cubic', help='For the geo2tiff step, define the resampling method (default: cubic). Options: near, bilinear, cubic, etc.')
    key_params.add_argument('--blocksize', type=int, default='512', help='For the geo2tiff step, define the blocksize (default: 512)')
    #key_params.add_argument('--aws-dir', type=str, default=None, help='For the update-aws step, define the path to the AWS S3 directory')
    key_params.add_argument('--catalog-yaml', type=str, default=None, help='For the update-catalog step, define the catalog yaml file to update (default: catalog.yaml in the output directory specified by --out-prefix)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles')
    aux_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    aux_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binary')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default=None, help='The file name of the Makefile to generate (default: {out-prefix}.mk)')
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it (default: False)')

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
    # - geotiff 
    geotif_f = args.in_fig 

    # - mbtiles
    mbtile_f=f"{args.out_prefix}.pmtiles.mbtiles"
    mbtile_f_ann=f"{args.out_prefix}.pmtiles.{args.resample}.mbtiles"
    # - pmtiles
    pmtiles_f=f"{args.out_prefix}.pmtiles"
    catalog_f = f"{out_dir}/catalog.yaml" if args.catalog_yaml is None else args.catalog_yaml

    # start mm
    mm = minimake()

    # 0. Create a geotiff file from PNG or TIF file
    if args.georeference:
        georef_f =f"{args.out_prefix}.georef.tif"
        assert args.in_fig is not None and os.path.exists(args.in_fig), "Please provide a valid input figure file using --in-figure"
        if args.in_bounds is not None:
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
        
        cmds = cmd_separator([], f"Geo-referencing {args.in_fig} to {georef_f}")
        cmds.append(f"{args.gdal_translate} -of GTiff -a_srs {args.srs} -a_ullr {ulx} {uly} {lrx} {lry} {args.in_fig} {georef_f}")
        mm.add_target(georef_f, [args.in_fig], cmds)
        # update geotif_f
        geotif_f = georef_f

    # 1. Orientation
    if args.flip_vertical or args.flip_horizontal or args.rotate is not None:
        # 0) Extract dimensions
        dim_f=geotif_f.replace(".tif","") + ".dim.tsv"
        cmds=cmds_for_dimensions(geotif_f, dim_f)
        mm.add_target(dim_f, [geotif_f], cmds)
    
        # 1) Output FILENAME
        ort_suffix=[]
        if args.rotate is not None:
            ort_suffix.append(f"rot{args.rotate}")
        if args.flip_vertical :
            ort_suffix.append("vflip")
        if args.flip_horizontal:
            ort_suffix.append("hflip")

        ort_f = f"{args.out_prefix}-{"-".join(ort_suffix)}.tif"

        # 2) Order argument
        # if args.rotate == "90" and not args.flip_vertical and not args.flip_horizontal:
        #     order_arg = "2,1"
        # elif args.rotate == "180" and not args.flip_vertical and not args.flip_horizontal:
        #     order_arg = "-1,-2"  
        # elif args.rotate is None and args.flip_vertical and args.flip_horizontal:  
        #     order_arg = "-1,-2"  
        # elif args.rotate == "270" and not args.flip_vertical and not args.flip_horizontal:
        #     order_arg = "2,-1"
        # elif args.flip_horizontal and not args.flip_vertical and args.rotate is None:
        #     order_arg = "-1,2"
        # elif args.flip_vertical and not args.flip_horizontal and args.rotate is None:
        #     order_arg = "1,-2"
        # elif args.rotate == "90" and args.flip_horizontal and not args.flip_vertical:
        #     order_arg = "-2,1"
        # elif args.rotate == "90" and args.flip_vertical and not args.flip_horizontal:
        #     order_arg = "-2,-1"
        # else:
        #     raise ValueError("Invalid combination of rotation and flip options.")

        order_arg = rotation_flip_map.get(
            (args.rotate, args.flip_vertical, args.flip_horizontal)
        )

        if order_arg is None:
            raise ValueError("Invalid combination of rotation and flip options.")

        # 3) Output Dimensions
        if order_arg.startswith("1") or order_arg.startswith("-1"):
            out_dim="$WIDTH $HEIGHT"
        else:
            out_dim="$HEIGHT $WIDTH"
        
        # 4) cmds
        msg=" ".join([f"Orientate the geotif, including ",
                "vertical flip;" if args.flip_vertical else "",
                "horizontal flip;" if args.flip_horizontal else "",
                f"rotate {args.rotate} degree clockwise;" if args.rotate is not None else "",
                geotif_f])
        cmds = cmd_separator([], msg) 
        cmds.append(f"WIDTH=$(awk '/WIDTH/' {dim_f}|cut -f 2) && \\")
        cmds.append(f"HEIGHT=$(awk '/HEIGHT/' {dim_f}|cut -f 2) && \\")
        cmd = " ".join([
                "gdalwarp",
                f'"{geotif_f}"',  # Add quotes around file names to handle spaces
                f'"{ort_f}"',
                "-b 1 -b 2 -b 3" if not args.mono else "-b 1",
                "-ct", f"\"+proj=pipeline +step +proj=axisswap +order={order_arg}\"",
                "-overwrite",
                "-ts", f"{out_dim}"
            ])
        cmds.append(cmd)
        mm.add_target(ort_f, [geotif_f, dim_f], cmds)

        # update geotif_f to the oriented file
        geotif_f = ort_f

    # 2. Convert a geotiff file to mbtiles
    if args.geotif2mbtiles:
        cmds = cmd_separator([], f"Converting from geotif to mbtiles: {geotif_f}")
        cmd = " ".join([
            "gdal_translate", 
            "-b 1 -b 2 -b 3" if not args.mono else "-b 1",
            "-strict",
            "-co", "\"ZOOM_LEVEL_STRATEGY=UPPER\"",
            "-co", f"\"RESAMPLING={args.resample}\"",
            "-co", f"\"BLOCKSIZE={args.blocksize}\"",
            "-ot", "Byte",
            "-scale", 
            "-of", "mbtiles",
            "-a_srs", args.srs,
            geotif_f, 
            mbtile_f
        ])
        cmds.append(cmd)
        mm.add_target(mbtile_f, [geotif_f], cmds)
    
    # 3. Convert mbtiles to pmtiles
    if args.mbtiles2pmtiles:
        cmds = cmd_separator([], f"Resampling mbtiles and converting to pmtiles: {geotif_f}")
        cmds.append(f"cp {mbtile_f} {mbtile_f_ann}")
        cmds.append(f"'{args.gdaladdo}' {mbtile_f_ann} -r {args.resample} 2 4 8 16 32 64 128 256")
        cmds.append(f"'{args.pmtiles}' convert --force {mbtile_f_ann} {pmtiles_f}")
        mm.add_target(pmtiles_f, [mbtile_f], cmds)

    # 4. Update the catalog.yaml file with the new pmtiles and upload to AWS
    if args.update_catalog:
        cmds = cmd_separator([], f"Updating yaml: {geotif_f}")
        pmtiles_name=os.path.basename(pmtiles_f)
        pmtiles_dir=os.path.dirname(pmtiles_f)
        cmds.append(f"cartloader update_catalog_for_basemap --in-yaml {catalog_f} --basemap {args.basemap_key}:{pmtiles_name} --basemap-dir {pmtiles_dir} --overwrite")
        cmds.append(f"touch {pmtiles_f}.yaml.done")
        mm.add_target(f"{pmtiles_f}.yaml.done", [pmtiles_f, catalog_f], cmds)

    # 5. Upload new PMtiles to AWS
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
        os.system(f"make -f {make_f} -n")
        print(f"To execute the pipeline, run the following command:\nmake -f {make_f} -j {args.n_jobs}")
    else:
        result = subprocess.run(f"make -f {make_f} -j {args.n_jobs} {'-B' if args.restart else ''}", shell=True)
        if result.returncode != 0:
            print(f"Error in converting the figure ({args.in_fig}) to pmtiles ({pmtiles_f})")
            sys.exit(1)

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])