import sys, os, argparse, inspect, subprocess, json

from cartloader.image import configure_color_mode
from cartloader.utils.utils import scheck_actions, write_dict_to_file, execute_makefile, cmd_separator, scheck_app, create_custom_logger, read_minmax, flexopen
from cartloader.utils.orient_helper import update_orient
import tifffile
import numpy as np


def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", description="""
                                    Convert Xenium H&E image to PMTiles format based on the alignment information using affine transformation.
                                     """)

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument('--dry-run', action='store_true', default=False, help='Dry run. Generate only the Makefile without running it')
    run_params.add_argument('--restart', action='store_true', default=False, help='Restart the run. Ignore all intermediate files and start from the beginning')
    run_params.add_argument('--n-jobs', type=int, default=1, help='Number of jobs (processes) to run in parallel')
    run_params.add_argument('--makefn', type=str, default=None, help='The file name of the Makefile to generate (default: {out-prefix}.png2pmtiles.mk)')

    inout_params = parser.add_argument_group("Input/Output Parameters")
    inout_params.add_argument('--tif', type=str, required=True, help='Path to the input H&E OME-TIFF file')
    inout_params.add_argument('--csv', type=str, required=True, help='Path to the Xenium alignment CSV file')
    inout_params.add_argument('--out-prefix', type=str, help='Prefix for naming the output files.')
    inout_params.add_argument('--px-size-xy', type=float, default=None, help='Override pixel size in um for both X and Y dimensions. If provided, skip reading OME-TIFF metadata.')

    key_params = parser.add_argument_group("Key parameters")
    key_params.add_argument('--threads', type=int, default=12, help='Number of threads to use for processing (default: 12)')
    key_params.add_argument('--remove-intermediate-files', action='store_true', default=False, help='If set, remove intermediate files (e.g., .mbtiles) after generating the final output.')
    
    env_params = parser.add_argument_group("Env Parameters", "Environment parameters, e.g., tools.")
    env_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='Path to pmtiles binary from go-pmtiles (default: pmtiles)')
    env_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary (default: gdal_translate)')
    env_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='Path to gdaladdo binar (default: gdaladdo)')
    env_params.add_argument('--gdalinfo', type=str, default=f"gdalinfo", help='Path to gdalinfo binary (default: gdalinfo)')
    env_params.add_argument('--vips', type=str, default=f"vips", help='Path to vips binary (default: vips)')

    aux_params = parser.add_argument_group("Auxiliary Parameters")    
    aux_params.add_argument('--log', action='store_true', default=False)
    aux_params.add_argument('--log-suffix', type=str, default=".log")

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def image_xenium_hne2pmtiles(_args):

    args=parse_arguments(_args)

    # dirs/files
    # - output dir
    out_dir = os.path.dirname(args.out_prefix)
    if len(out_dir) > 0:
        print(f"Creating output directory: {out_dir}")
        os.makedirs(out_dir, exist_ok=True)
    
    # - input image
    assert os.path.exists(args.tif), f"File not found: {args.tif} (--tif)"
    assert os.path.exists(args.csv), f"File not found: {args.csv} (--csv)"

    logger = create_custom_logger(__name__, args.out_prefix + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    ## Read the CSV file
    trans_mat = np.zeros((3,3), dtype=np.float64)
    with flexopen(args.csv, 'rt') as rf:
        toks = rf.readline().split(',')
        trans_mat[0,0] = float(toks[0])
        trans_mat[0,1] = float(toks[1])
        trans_mat[0,2] = float(toks[2])

        toks = rf.readline().split(',')
        trans_mat[1,0] = float(toks[0])
        trans_mat[1,1] = float(toks[1])
        trans_mat[1,2] = float(toks[2])

        toks = rf.readline().split(',')
        trans_mat[2,0] = float(toks[0])
        trans_mat[2,1] = float(toks[1])
        trans_mat[2,2] = float(toks[2])

    logger.info("Transformation Matrix: \n" + str(trans_mat))
    if trans_mat[1,0] > 0 or trans_mat[0,1] < 0:
        raise ValueError("The transformation matrix seems incorrect. Please check the CSV file.")
    scale = np.sqrt(np.linalg.det(trans_mat))

    if args.px_size_xy is not None:
        px_size_xy = args.px_size_xy
        logger.info(f"Using overridden pixel size: {px_size_xy} um")
    else:
        logger.info("Reading OME-TIFF metadata for pixel size...")
        # [CSV processing and metadata extraction remain the same]
        with tifffile.TiffFile(args.tif, _multifile=False) as tif:
            logger.info(f"Loaded OME-TIFF file {args.tif}")
            
            meta = tifffile.xml2dict(tif.ome_metadata) ## extract metadata
            meta = meta['OME']['Image']['Pixels']
            logger.info("Metadata: \n" + json.dumps(meta, indent=2))            
            px_size_x = float(meta['PhysicalSizeX'])
            px_size_y = float(meta['PhysicalSizeY'])
            px_size_unit = meta.get('PhysicalSizeXUnit', 'um')
            if px_size_unit != 'um' and px_size_unit != 'µm':
                raise ValueError(f"Physical size unit is not in um: {px_size_unit}")
            offset_um_x = float(meta.get('OffsetX', 0))
            offset_um_y = float(meta.get('OffsetY', 0))

            logger.info("Physical Size X: {:.4f} um".format(px_size_x))
            logger.info("Physical Size Y: {:.4f} um".format(px_size_y))
            logger.info("Offset X: {:.4f} um".format(offset_um_x))
            logger.info("Offset Y: {:.4f} um".format(offset_um_y))

            ## Generate affine transformation matrix
            px_size_xy = (px_size_x + px_size_y) / 2.0

    res_ratio = scale / px_size_xy
    affine_mat = trans_mat * res_ratio

    vec_x = affine_mat[0:2, 0]
    vec_y = affine_mat[0:2, 1]

    affine_mat[0, 2] -= 0.5 * (vec_x[0] + vec_y[0])
    affine_mat[1, 2] -= 0.5 * (vec_x[1] + vec_y[1])
    #affine_mat = trans_mat / scale * px_size_xy
    
    logger.info("Affine Transformation Matrix (in um): \n" + str(affine_mat))

    cmd = f"GDAL_NUM_THREADS={args.threads} {args.gdal_translate} -of GTiff -a_srs EPSG:3857 -a_gt {affine_mat[0,2]} {affine_mat[0,0]} {affine_mat[0,1]} {affine_mat[1,2]} {affine_mat[1,0]} {affine_mat[1,1]} '{args.tif}' '{args.out_prefix}.georef.tif'"
    logger.info("Running gdal_translate command:\n" + cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("gdal_translate failed with error:\n" + result.stderr)
        raise RuntimeError("gdal_translate command failed.")
        # fallback 
        # logger.info("Trying vips as a fallback...")
        # cmd = f"{args.vips} copy '{args.tif}' '{args.out_prefix}.png'"
        # logger.info("Running vips command:\n" + cmd)
        # result2 = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        # if result2.returncode != 0:
        #     logger.error("vips failed with error:\n" + result2.stderr)
        #     raise RuntimeError("Both gdal_translate and vips commands failed.")
        # cmd = f"GDAL_NUM_THREADS={args.threads} {args.gdal_translate} -of GTiff -a_srs EPSG:3857 -a_gt {affine_mat[0,2]} {affine_mat[0,0]} {affine_mat[0,1]} {affine_mat[1,2]} {affine_mat[1,0]} {affine_mat[1,1]} '{args.out_prefix}.png' '{args.out_prefix}.georef.tif'"
        # logger.info("Running gdal_translate command:\n" + cmd)
        # result2 = subprocess.run(cmd, shell=True, capture_output=True, text=True)        
        # if result2.returncode != 0:
        #     logger.error("gdal_translate failed with error:\n" + result2.stderr)
        #     raise RuntimeError("gdal_translate command failed even after vips conversion.")
    
    cmd = f"GDAL_NUM_THREADS={args.threads} gdal_translate -ot Byte -of mbtiles -b 1 -b 2 -b 3 -strict -co 'ZOOM_LEVEL_STRATEGY=UPPER' -co 'RESAMPLING=cubic' -co 'BLOCKSIZE=512' -co 'QUALITY=100' -a_srs EPSG:3857 {args.out_prefix}.georef.tif {args.out_prefix}.mbtiles"
    logger.info("Running gdal_translate command:\n" + cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("gdal_translate failed with error:\n" + result.stderr)
        raise RuntimeError("gdal_translate command failed.")
    
    cmd = f"GDAL_NUM_THREADS={args.threads} {args.gdaladdo} '{args.out_prefix}.mbtiles' -r cubic 2 4 8 16 32 64 128 256 512 1024"
    logger.info("Running gdal_translate command:\n" + cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("gdal_translate failed with error:\n" + result.stderr)
        raise RuntimeError("gdaladdo command failed.")

    cmd = f"{args.pmtiles} convert --force '{args.out_prefix}.mbtiles' '{args.out_prefix}.pmtiles'"
    logger.info("Running pmtiles convert command:\n" + cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("pmtiles convert failed with error:\n" + result.stderr)
        raise RuntimeError("pmtiles convert command failed.")
    
    if args.remove_intermediate_files:
        os.remove(f"{args.out_prefix}.mbtiles")
        os.remove(f"{args.out_prefix}.georef.tif")
        logger.info("Removed intermediate files.")
    logger.info("Analysis Finished")


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
