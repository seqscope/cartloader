import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast
import pandas as pd

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger


def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

    parser = argparse.ArgumentParser(prog=f"cartloader run_tsv2pmtiles", description="Split and convert transcripts TSV file into pmtiles")

    cmd_params = parser.add_argument_group("Commands", "Commands to run together")
    cmd_params.add_argument('--all', action='store_true', default=False, help='Run all commands (geo2tiff, mbtiles2pmtiles, resample, mbtiles2pmtiles)')
    cmd_params.add_argument('--georeference', action='store_true', default=False, help='Create a geotiff file from PNG or TIF file. If enabled, the user must provide georeferenced bounds using --in-tsv or --in-bounds.')
    cmd_params.add_argument('--geotiff2mbtiles', action='store_true', default=False, help='Convert a geotiff file to mbtiles')
    cmd_params.add_argument('--mbtiles2pmtiles', action='store_true', default=False, help='Convert mbtiles to pmtiles')
    cmd_params.add_argument('--upload-aws', action='store_true', default=False, help='Upload the new pmtiles and updated catalog.yaml file to the AWS S3 bucket')

    inout_params = parser.add_argument_group("Input/Output Parameters", "Define the input file according to the user's needs.")
    inout_params.add_argument('--in-fig', type=str, help='The input figure file (PNG or TIF) to be converted to pmTiles')
    inout_params.add_argument('--in-tsv', type=str, default=None, help='If --georeference is required, use the *.pixel.sorted.tsv.gz from run_ficture to provide georeferenced bounds.')
    inout_params.add_argument('--in-bounds', type=str, default=None, help='If --georeference is required, provide the bounds in the format of "<ulx>,<uly>,<lrx>,<lry>", which represents upper-left X, upper-left Y, lower-right X, lower-right Y.')
    inout_params.add_argument('--out-prefix', required= True, type=str, help='The output prefix. New directory will be created if needed')

    key_params = parser.add_argument_group("Key parameters")
    key_params.add_argument('--srs', type=str, default='EPSG:3857', help='For the georeference and geo2tiff steps, define the spatial reference system (default: EPSG:3857)')
    key_params.add_argument('--resample-method', type=str, default='nearest', help='For the geo2tiff step, define the resampling method (default: nearest)')
    key_params.add_argument('--blocksize', type=int, default='512', help='For the geo2tiff step, define the blocksize (default: 512)')
    key_params.add_argument('--aws-bucket', type=str, default=None, help='For the update-catalog-aws step, define the path to the AWS S3 bucket')
    
    key_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    key_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Auxiliary parameters (using default is recommended)")
    aux_params.add_argument('--tippecanoe', type=str, default=f"{repo_dir}/submodules/tippecanoe/tippecanoe", help='Path to tippecanoe binary')
    aux_params.add_argument('--pmtiles', type=str, default=f"{repo_dir}/submodules/go-pmtiles/pmtiles", help='Path to pmtiles binary from go-pmtiles')

    run_params = parser.add_argument_group("Run Options", "Run options for FICTURE commands")
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    #run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def run_fig2pmtiles(_args):

    args=parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_prefix + "_tsv2pmtiles" + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    if args.all:
        args.georeference = True
        args.geotiff2mbtiles = True
        args.mbtiles2pmtiles = True
        args.upload_aws = True

    # files
    out_dir = os.path.dirname(args.out_prefix)
    os.makedirs(out_dir, exist_ok=True)
    # - geotiff
    geotiff_f = args.in_fig if not args.georeference else f"{args.out_prefix}.pmtiles.tif"
    # - mbtiles
    mbtile_f=f"{args.out_prefix}.pmtiles.mbtiles"
    mbtile_f_ann=f"{args.out_prefix}.pmtiles.{args.resample_method}.mbtiles"
    # - pmtiles
    pmtiles_f=f"{args.out_prefix}.pmtiles"

    # start mm
    mm = minimake()

    # 1. Create a geotiff file from PNG or TIF file
    if args.georeference:
        assert args.in_fig is not None and os.path.exists(args.in_fig), "Please provide a valid input figure file using --in-figure"
        if args.in_bounds is not None:
            ulx,uly,lrx,lry = args.in_bounds.split(",")
        elif args.in_tsv is not None:
            with gzip.open(args.in_tsv, 'rt') as f:
                for i in range(3):
                    line = f.readline()
            ann2val = {x.split("=")[0]:x.split("=")[1] for x in line.strip().replace("##", "").split(";")}
            ulx = float(ann2val["OFFSET_X"])
            uly = float(ann2val["OFFSET_Y"])
            lrx = float(ann2val["SIZE_X"])+1+ulx
            lrx = float(ann2val["SIZE_Y"])+1+uly
        else:
            raise ValueError("Please provide either --in-bounds or --in-tsv to georeference the figure")
        
        cmds = cmd_separator([], f"Geo-referencing {args.in_fig} to {geotiff_f}")
        cmds.append(f"gdal_translate -of GTiff -a_srs {args.srs} -a_ullr {ulx} {uly} {lrx} {lry} {args.in_fig} {geotiff_f}")
        mm.add_target(geotiff_f, [args.in_fig], cmds)

    # 2. Convert a geotiff file to mbtiles
    if args.geotiff2mbtiles:
        cmds = cmd_separator([], f"Converting a geotif to mbtiles")
        cmd = " ".join([
            "gdal_translate", 
            "-b", "1", 
            "-b", "2",
            "-b", "3",
            "-strict",
            "-co", "\"ZOOM_LEVEL_STRATEGY=UPPER\"",
            "-co", "\"RESAMPLING={args.resample_method}\"",
            "-co", "\"BLOCKSIZE={args.blocksize}\"",
            "-ot", "Byte",
            "-scale", "0", "255",
            "-of", "mbtiles",
            "-a_srs", args.srs,
            "geotiff_f", mbtile_f
        ])
        mm.append(cmd)
        #cmds.append(f"gdal_translate -b 1 -b 2 -b 3 -strict -co \"ZOOM_LEVEL_STRATEGY=UPPER\" -co \"RESAMPLING={args.resample_method}\" -co \"BLOCKSIZE={args.blocksize}\" -ot Byte -scale 0 255 -of mbtiles -a_srs {args.srs} {geotiff_f} {mbtile_f}")
        mm.add_target(mbtile_f, [geotiff_f], cmds)
    
    # 3. Convert mbtiles to pmtiles
    if args.mbtiles2pmtiles:
        cmds = cmd_separator([], f"Resampling mbtiles and converting to pmtiles")
        cmds.append(f"cp {mbtile_f} {mbtile_f_ann}")
        cmds.append(f"gdaladdo {mbtile_f_ann} -r {args.resample_method} 2 4 8 16 32 64 128 256")
        cmds.append(f"{args.tippecanoe} convert --force {mbtile_f_ann} {pmtiles_f}")
        cmds.append(f"rm {mbtile_f_ann}")
        mm.add_target(pmtiles_f, [mbtile_f], cmds)

    # 4. Update to catalog and upload to AWS
    if args.update_aws:
        cmds = cmd_separator([], f"Updating catalog.yaml")
        cmds.append(f"cartloader update_yaml_for_basemap --catalog {args.catalog} --pmtiles {pmtiles_f}")
        cmds.append(f"aws s3 cp {pmtiles_f} {args.aws_bucket}")
        cmds.append(f"aws s3 cp {args.catalog} {args.aws_bucket}/catalog.yaml")
        # touch a done file
        cmds.append(f"touch {pmtiles_f}.done")
        mm.add_target(f"{out_dir}/{pmtiles_f}.done", [pmtiles_f], cmds)
