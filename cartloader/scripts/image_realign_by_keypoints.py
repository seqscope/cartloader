import sys, os, argparse, logging,  inspect, json, subprocess
import pandas as pd
import shutil
from cartloader.utils.utils import flexopen, read_minmax, create_custom_logger
from cartloader.utils.image_helper import check_north_up

# get the path of the cu

def parse_arguments(_args):
    """Parse command-line arguments."""
    repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""
                                    Realign a GeoTIFF image based on the list of key points.
                                     """)
    parser.add_argument('--in-tif', type=str, required=True,help='Input GeoTIFF file. Must be already georeferenced')
    parser.add_argument('--out', type=str, required=True, help='Output prefix')
    parser.add_argument('--keypoints', type=str, required=True, help='Input keypoints file [new_X] [new_Y] [old_X] [old_Y] in each line in EPSG:3857 coordinate')
    parser.add_argument('--gdaltransform', type=str, default=f"gdaltransform", help='Path to gdal_transform binary')
    parser.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='Path to gdal_translate binary')
    parser.add_argument('--gdalwarp', type=str, default=f"gdalwarp", help='Path to gdalwarp binary')
    parser.add_argument('--geotiff2pmtiles', type=str, default=f"{repo_dir}/submodules/geotiff2pmtiles/geotiff2pmtiles", help='Path to geotiff2pmtiles binary (default: geotiff2pmtiles)')
    parser.add_argument('--min-zoom', type=int, default=6, help='Minimum zoom level for PMTiles (default: 6)')
    parser.add_argument('--max-zoom', type=int, help='Maximum zoom level for PMTiles (default: auto)')
    parser.add_argument('--tile-format', type=str, default='png', choices=['png', 'webp'], help='Tile format for PMTiles (default: png)')
    parser.add_argument('--log', action='store_true', default=False, help='Write logs to a file under the output directory')
    parser.add_argument('--log-suffix', type=str, default=".log", help='Suffix for the log filename; final path is <out_dir>_cartload<suffix> (default: .log)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def image_realign_by_keypoints(_args):
    args = parse_arguments(_args)

    logger = create_custom_logger(__name__, args.out_prefix + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    out_dir = os.path.dirname(args.out)
    if out_dir != "":
        os.makedirs(out_dir, exist_ok=True)

    ## write trans.txt file first
    logger.info(f"Writing transformation keypoints to {args.out}.src.txt")
    new_xys = []
    with flexopen(args.keypoints, "r") as f:
        with flexopen(f"{args.out}.src.txt", "w") as wf:
            for line in f:
                newx, newy, oldx, oldy = line.strip().split()
                new_xys.append((newx, newy))
                wf.write(f"{oldx} {oldy}\n")

    ## run gdaltransform
    logger.info(f"Running gdaltransform to compute the transformation parameters")
    cmd = f"cat {args.out}.src.txt | {args.gdaltransform} -i {args.in_tif} > {args.out}.trans.txt"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    logger.info(f"running gdaltransform command: {cmd}")
    if result.returncode != 0:
        logger.error(f"Error running gdaltransform: {result.stderr}")
        raise RuntimeError(f"gdaltransform failed with error: {result.stderr}")
    logger.info(f"gdaltransform completed successfully. Transformation parameters written to {args.out}.trans.txt")

    ## read the transformed keypoints
    src_xys = []
    with flexopen(f"{args.out}.trans.txt", "r") as f:
        for line in f:
            srcx, srcy, z = line.strip().split()
            src_xys.append((srcx, srcy))

    ## create VRT file
    logger.info(f"Creating VRT file for alignment")
    cmd = f"{args.gdal_translate} -of VRT "
    for i in range(len(new_xys)):
        cmd += f" -gcp {src_xys[i][0]} {src_xys[i][1]} {new_xys[i][0]} {new_xys[i][1]}"
    cmd += f" {args.in_tif} {args.out}.vrt"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    logger.info(f"running gdal_translate command for VRT creation: {cmd}")
    if result.returncode != 0:
        logger.error(f"Error running gdal_translate for VRT creation: {result.stderr}")
        raise RuntimeError(f"gdal_translate for VRT creation failed with error: {result.stderr}")
    logger.info(f"VRT file created successfully at {args.out}.vrt")

    ## run gdalwarp to realign the image
    logger.info(f"Running gdalwarp to realign the image")
    cmd = f"{args.gdalwarp} -tps -r bilinear -t_srs EPSG:3857 {args.out}.vrt {args.out}.tif"
    logger.info(f"running gdalwarp command for image realignment: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Error running gdalwarp for image realignment: {result.stderr}")
        raise RuntimeError(f"gdalwarp for image realignment failed with error: {result.stderr}")
    logger.info(f"Image realignment completed successfully. Realigned image written to {args.out}.tif")


    ## run geotiff2pmtiles to generate pmtiles
    logger.info(f"Running geotiff2pmtiles to generate PMTiles")
    cmd = f"{args.geotiff2pmtiles} --format {args.tile_format} --min-zoom {args.min_zoom} "
    if args.max_zoom is not None:
        cmd += f"--maxzoom {args.max_zoom} "
    cmd += f"{args.out}.tif {args.out}.pmtiles"
    logger.info(f"running geotiff2pmtiles command for PMTiles generation: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Error running geotiff2pmtiles for PMTiles generation: {result.stderr}")
        raise RuntimeError(f"geotiff2pmtiles for PMTiles generation failed with error: {result.stderr}")
    logger.info(f"PMTiles generation completed successfully. PMTiles written to {args.out}.pmtiles")


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
