import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, re
import pandas as pd
import numpy as np
import tifffile, json
from PIL import Image

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, read_minmax, hex_to_rgb

def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(prog=f"cartloader transform_aligned_histology", description="Convert Aligned Histology from 10x Xenium or Vizgen MERSCOPE")

    inout_params = parser.add_argument_group("Input/Output Parameters", "Define the input file according to the user's needs.")
    inout_params.add_argument('--tif', type=str, required=True, help='Input OME-TIFF file')
    inout_params.add_argument('--csv', type=str, help='(Vizgen only) CSV file containing conversion parameters (typically micron_to_mosaic_pixel_transform.csv)')
    inout_params.add_argument('--out-prefix', type=str, required=True, help='Output file prefix')
    inout_params.add_argument("--page", type=int, help='For 3D (X/Y/Z) OME file, specify the z-value to extract the image from')
    inout_params.add_argument("--level", type=int, default=0, help='Specify the level to extract from the OME-TIFF file')
    inout_params.add_argument("--series", type=int, default=0, help='Specify the index of series to extract from the OME-TIFF file')
    inout_params.add_argument("--upper-thres-quantile", type=float, help='Quantile-based capped value for rescaling the image. Cannot be used with --upper-thres-intensity')
    inout_params.add_argument("--upper-thres-intensity", type=float, default=255, help='Intensity-based capped value for rescaling the image. Cannot be used with --upper-thres-quantile')
    inout_params.add_argument("--lower-thres-quantile", type=float, help='Quantile-based floored value for rescaling the image. Cannot be used with --lower-thres-intensity')
    inout_params.add_argument("--lower-thres-intensity", type=float, default=0, help='Intensity-based floored value for rescaling the image. Cannot be used with --lower-thres-quantile')
    inout_params.add_argument("--colorize", type=str, help='Colorize the black-and-white image using a specific RGB code as a max value (does not work with RGB images)')
    inout_params.add_argument('--flip-horizontal', action='store_true', default=False, help='Create PMTiles in addition to png')
    inout_params.add_argument('--flip-vertical', action='store_true', default=False, help='Create PMTiles in addition to png')
    inout_params.add_argument('--rotate-clockwise', action='store_true', default=False, help='Rotate by 90 degree by clockwise direction (before flip)')
    inout_params.add_argument('--rotate-counter', action='store_true', default=False, help='Rotate by 90 degree by counterclockwise direction (before flip)')

    aux_params = parser.add_argument_group("Auxiliary Parameters", "Additional parameters for the script")    
    aux_params.add_argument('--skip-pmtiles', action='store_true', default=False, help='Skip creating PMTiles in addition to png')
    aux_params.add_argument('--write-mono', action='store_true', default=False, help='Write the is_mono flag in a separate file')
    aux_params.add_argument('--write-bounds', action='store_true', default=False, help='Write the bounds of <ulx>,<uly>,<lrx>,<lry> in a separate file')
    aux_params.add_argument('--log', action='store_true', default=False, help='Write log to file')
    aux_params.add_argument('--log-suffix', type=str, default=".log", help='The suffix for the log file (appended to the output directory). Default: .log')
    aux_params.add_argument('--pmtiles', type=str, default=f"pmtiles", help='If not --skip-pmtile, path to pmtiles binary from go-pmtiles')
    aux_params.add_argument('--gdal_translate', type=str, default=f"gdal_translate", help='If not --skip-pmtile, path to gdal_translate binary')
    aux_params.add_argument('--gdaladdo', type=str, default=f"gdaladdo", help='If not --skip-pmtile, path to gdaladdo binary')
    # aux_params.add_argument('--keep-intermediate-files', action='store_true', default=False, help='Keep intermediate files')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def transform_aligned_histology(_args):
    args=parse_arguments(_args)
        
    ## ~/data/ficture/xenium/ffpe_human_breast_cancer_rep1/out/outs/morphology.ome.tif

    logger = create_custom_logger(__name__, args.out_prefix + args.log_suffix if args.log else None)
    logger.info("Analysis Started")
    
    # replace is_ome = True with check_ome_tiff ??
    # def check_ome_tiff(tiff_path):
    #     assert os.path.exists(tiff_path), f"The file '{tiff_path}' does not exist."
    #     with tifffile.TiffFile(tiff_path) as tif:
    #         if tif.is_ome:
    #             return True
    #         else:
    #             return False
    is_ome = True
    px_per_um_x = None
    px_per_um_y = None
    offset_px_x = None
    offset_px_y = None
    if args.csv is not None:
        is_ome = False
        logger.info(f"Reading the CSV file {args.csv}")
        with open(args.csv, 'r') as f:
            lines = f.readlines()
            tok0s = lines[0].strip().replace(',', ' ').split(' ')
            tok1s = lines[1].strip().replace(',', ' ').split(' ')
            px_per_um_x = float(tok0s[0])
            px_per_um_y = float(tok1s[1])
            offset_px_x = float(tok0s[2])
            offset_px_y = float(tok1s[2])
            if float(tok0s[1]) != 0.0 or float(tok1s[0]) != 0.0:
                raise ValueError(f"The CSV file {args.csv} contains non-zero values at off-diagonal elements")
            
    ## output basic information, such as number of pages
    with tifffile.TiffFile(args.tif, _multifile=False) as tif:
        logger.info(f"Successfully loaded the OME-TIFF file {args.tif}")
        n_pages = len(tif.pages)
        n_series = len(tif.series)
        # ?? Why using tif.series[0] throughout the script while args.series exists?
        # ?? replace it with `n_levels = len(tif.series[args.series].levels)` given `page = tif.series[args.series].levels[args.level].pages[args.page]`
        n_levels = len(tif.series[0].levels)
        
        # ??
        # args.page = 0
        if args.page is None:
            if n_pages > 1:
                logger.error("Multiple pages detected. Please specify the page number to extract the image from")
                sys.exit(1)
            elif n_pages == 1:
                args.page = 0
            else:
                raise ValueError("In --transform, no pages detected in the OME-TIFF file")

        page = tif.series[args.series].levels[args.level].pages[args.page]

## This is a test code for handle the tiled TIFF      
        if page.is_tiled:            
            print(f"This is a tiled TIFF with chunk size {page.chunks}")
            print(f"chunks = {page.chunks}")
            print(f"tilelength = {page.tilelength}")
            print(f"tilewidth = {page.tilewidth}")
            print(f"tiledepth = {page.tiledepth}")
            print(f"samplesperpixel = {page.samplesperpixel}")
            print(f"imagewidth = {page.imagewidth}")
            print(f"imagelength = {page.imagelength}")
            print(f"imagedepth = {page.imagedepth}")
            print(f"rowsperstrip = {page.rowsperstrip}")
            print(f"sampleformat = {page.sampleformat}")
            print(f"compression = {page.compression}")
            print(f"photometric = {page.photometric}")
            segments = page.segments()
            for i, segment in enumerate(segments):
                (data, offset, bytecount) = segment
                print(f"{i}, {data.shape}, {offset}, {bytecount}, {np.squeeze(data).shape}")
                if i > 10:
                    break
        else:
            # It's not tiled, so check how many strips there are
            n_strips = len(page.strip_offsets) if page.strip_offsets is not None else 0
            if n_strips > 1:
                print(f"This is a striped TIFF with {n_strips} strips.")
            else:
                print("This is a single-strip (scanline-based) TIFF.")
        
        sys.exit(1)
            
        if len(page.shape) == 3:
            if page.shape[2] != 3:
                logger.error("The colored image is not in RGB format")
                sys.exit(1)
            is_mono = False
        else:
            is_mono = True
            
        if args.colorize is not None and is_mono is False:
            logger.error("Cannot colorize the RGB-colored image")
            sys.exit(1)

        if is_ome:
            meta = tifffile.xml2dict(tif.ome_metadata) ## extract metadata
            meta = meta['OME']['Image']['Pixels']
            logger.info("Metadata: \n" + json.dumps(meta, indent=2))            
            px_size_x = meta['PhysicalSizeX']
            px_size_y = meta['PhysicalSizeY']
            px_size_unit = meta.get('PhysicalSizeXUnit', 'um')
            if px_size_unit != 'um' and px_size_unit != 'µm':
                raise ValueError(f"Physical size unit is not in um: {px_size_unit}")
            offset_um_x = meta.get('OffsetX', 0)
            offset_um_y = meta.get('OffsetY', 0)

            level_0_shape = tif.series[0].levels[0].pages[args.page].shape
            current_page_shape = page.shape
            if level_0_shape != current_page_shape:
                scale_factor_x = level_0_shape[1] / current_page_shape[1]
                scale_factor_y = level_0_shape[0] / current_page_shape[0]
                px_size_x = px_size_x * scale_factor_x
                px_size_y = px_size_y * scale_factor_y
                logger.info(f"Rescaling the pixel size level {args.level} by ({scale_factor_x},{scale_factor_y})...")
        else:
            px_size_x = 1/px_per_um_x
            px_size_y = 1/px_per_um_y
            offset_um_x = 0 - offset_px_x / px_per_um_x
            offset_um_y = 0 - offset_px_y / px_per_um_y
            
        ul = [offset_um_x, offset_um_y]
        lr = [offset_um_x + px_size_x * page.shape[1], offset_um_y + px_size_y * page.shape[0]]

        logger.info(f"OME-TIFF mode: {is_ome}")        
        logger.info(f"Number of pages: {n_pages}")
        logger.info(f"Number of series: {n_series}")
        logger.info(f"Number of levels: {n_levels}")
        logger.info(f"Pixel size in um: ({px_size_x},{px_size_y})")
        logger.info(f"Offset in um: ({offset_um_x},{offset_um_y})")
        for i in range(n_levels):
            logger.info(f"Level {i} Dimension: {tif.series[0].levels[i].shape}")
            #logger.info(type(tif.series[0].levels[i]))
                                        
        image = page.asarray()
        logger.info(f"Extracted image shape: {image.shape}")
        #logger.info(f"Extracted image type: {type(image.shape)}")
        
        if args.upper_thres_quantile is not None or args.lower_thres_quantile is not None:
            upper_thres_quantile = args.upper_thres_quantile if args.upper_thres_quantile is not None else 1
            lower_thres_quantile = args.lower_thres_quantile if args.lower_thres_quantile is not None else 0
            quantile_values = np.quantile(image, [0, 0.25, 0.5, 0.75, 1, upper_thres_quantile, lower_thres_quantile])
            logger.info(f"Minimum pixel intensity: {quantile_values[0]}")
            logger.info(f"Maximum pixel intensity: {quantile_values[4]}")
            logger.info(f"Median pixel intensity: {quantile_values[2]}")
            logger.info(f"IQR of pixel intensity: [{quantile_values[1]}, {quantile_values[3]}]")
            logger.info(f"{args.upper_thres_quantile} quantile threshold to cap the intensities: {quantile_values[5]}")
            logger.info(f"{args.lower_thres_quantile} quantile threshold to cap the intensities: {quantile_values[6]}")
            if args.upper_thres_quantile is not None:
                args.upper_thres_intensity = quantile_values[5]
            if args.lower_thres_quantile is not None:
                args.lower_thres_intensity = quantile_values[6]
        
        logger.info(f"Rescaling the image into 8 bits using [{args.lower_thres_intensity}, {args.upper_thres_intensity}] as clipping thresholds")
        image = ((np.clip(image, a_min=args.lower_thres_intensity, a_max=args.upper_thres_intensity) - args.lower_thres_intensity) / (args.upper_thres_intensity - args.lower_thres_intensity))
        if args.colorize is not None:
            r, g, b = hex_to_rgb(args.colorize)            
            image = np.stack([image * r, image * g, image * b], axis=-1).astype(np.uint8)
            is_mono = False
        else:
            image = (image * 255.0).astype(np.uint8)
        # #image[image > args.thres_intensity] = args.thres_intensity
        # if args.thres_intensity == 255:
        #     image = image.astype(np.uint8)            
        # else:
        #     image = (image * 255 / args.thres_intensity).astype(np.uint8)
        # #image = image.astype(np.uint8)
                
        if args.rotate_clockwise:
            logger.info(f"Performing 90-degree clockwise rotation in memory...")
            image = np.rot90(image, k=1, axes=(0, 1))
            (ul[0], ul[1], lr[0], lr[1]) = (ul[1], ul[0], lr[1], lr[0])
        if args.rotate_counter:
            logger.info(f"Performing 90-degree counterclockwise rotation in memory...")
            image = np.rot90(image, k=-1, axes=(0, 1))
            (ul[0], ul[1], lr[0], lr[1]) = (ul[1], ul[0], lr[1], lr[0])
        if args.flip_horizontal:
            logger.info(f"Performing horizontal flip in memory...")
            image = np.fliplr(image)                
        if args.flip_vertical:
            logger.info(f"Performing vertical flip in memory...")
            image = np.flipud(image)                

        logger.info(f"Converting the image to PNG using PIL..")
        #print(image)            
        image = Image.fromarray(image)
        
        logger.info(f"Saving the image...")
        image.save(f"{args.out_prefix}.png") 

        if args.write_bounds:
            if ul[0] < 0:
                ul0 = f"\\{ul[0]}"
            else:
                ul0 = f"{ul[0]}"
            with open(f"{args.out_prefix}.bounds.csv", 'w') as f:
                f.write(f"{ul[0]},{ul[1]},{lr[0]},{lr[1]}\n")

        if not args.skip_pmtiles:
            logger.info(f"Creating PMTiles with run_fig2pmtiles...")
            # if args.flip_horizontal:
            #     (ul[0], lr[0]) = (lr[0], ul[0])
            # if args.flip_vertical:
            #     (ul[1], lr[1]) = (lr[1], ul[1])
            if ul[0] < 0:
                ul0 = f"\\{ul[0]}"
            else:
                ul0 = f"{ul[0]}"
            cmd = f"cartloader run_fig2pmtiles --in-bounds '{ul0},{ul[1]},{lr[0]},{lr[1]}' --in-fig {args.out_prefix}.png --out-prefix {args.out_prefix} --geotif2mbtiles --mbtiles2pmtiles --georeference --pmtiles {args.pmtiles} --gdal_translate {args.gdal_translate} --gdaladdo {args.gdaladdo}"

            if is_mono:
                cmd += " --mono"
            print(cmd)
            result = subprocess.run(cmd, shell=True)
            if result.returncode != 0:
                print(f"Error in converting the PNG to pmtiles")
                sys.exit(1)
        
    logger.info("Analysis Completed")
    
if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])