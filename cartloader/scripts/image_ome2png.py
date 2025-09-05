import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, re, inspect
import pandas as pd
import numpy as np
import tifffile, json
from PIL import Image
import psutil
from typing import Tuple
from collections import Counter
import gc

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, read_minmax, hex_to_rgb

def parse_arguments(_args):
    # [Previous argument parsing code remains the same]
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",  description="Convert Aligned Histology from 10x Xenium or Vizgen MERSCOPE to a flat image")

    # Add existing arguments
    inout_params = parser.add_argument_group("Input/Output Parameters")
    inout_params.add_argument('--tif', type=str, required=True, help='Path to the input OME-TIFF file')
    inout_params.add_argument('--csv', type=str, help='(Vizgen only) Path to the input CSV file containing conversion parameters')
    inout_params.add_argument('--out-prefix', type=str, required=True, help='Output file prefix')
    inout_params.add_argument("--page", type=int, help='For 3D (X/Y/Z) OME file, specify the z-value')
    inout_params.add_argument("--level", type=int, default=0, help='Specify the level to extract (default: 0)')
    inout_params.add_argument("--series", type=int, default=0, help='Specify the index of series (default: 0)')
    inout_params.add_argument("--upper-thres-quantile", type=float, help='Quantile-based capped value')
    inout_params.add_argument("--upper-thres-intensity", type=float, default=255, help='Intensity-based capped value (default: 255)')
    inout_params.add_argument("--lower-thres-quantile", type=float, help='Quantile-based floored value')
    inout_params.add_argument("--lower-thres-intensity", type=float, default=0, help='Intensity-based floored value (default: 0)')
    inout_params.add_argument("--transparent-below", type=int, default=0, help='Set pixels below this value to transparent (0-255; default:0)')
    inout_params.add_argument("--colorize", type=str, help='Colorize using RGB code as max value')
    inout_params.add_argument('--flip-horizontal', action='store_true', default=False)
    inout_params.add_argument('--flip-vertical', action='store_true', default=False)
    inout_params.add_argument('--rotate-clockwise', action='store_true', default=False)
    inout_params.add_argument('--rotate-counter', action='store_true', default=False)
    inout_params.add_argument('--high-memory', action='store_true', default=False)
    inout_params.add_argument('--write-color-mode', action='store_true', default=False,  help='Save the color mode into a file (<out_prefix>.color.csv). This argument is specifically designed for "cartloader import_image"')

    # memory_params = parser.add_argument_group("Memory Management")
    # memory_params.add_argument('--max-memory-gb', type=float, default=4.0,
    #                          help='Maximum memory usage in gigabytes')

    aux_params = parser.add_argument_group("Auxiliary Parameters")    
    aux_params.add_argument('--log', action='store_true', default=False)
    aux_params.add_argument('--log-suffix', type=str, default=".log")

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(_args)

    # Validate transparent range
    if args.transparent_below is not None:
        if not (0 <= args.transparent_below <= 255):
            parser.error("--transparent-below must be between 0 and 255 inclusive")

    return args

def get_memory_usage():
    """Get current memory usage in GB"""
    process = psutil.Process()
    return process.memory_info().rss / (1024 * 1024 * 1024)

def process_image_chunk(chunk: np.ndarray, args, r=None, g=None, b=None) -> np.ndarray:
    """Process a chunk of the image with the specified transformations"""
    # Clip and normalize
    chunk = ((np.clip(chunk, a_min=args.lower_thres_intensity, a_max=args.upper_thres_intensity) 
             - args.lower_thres_intensity) / (args.upper_thres_intensity - args.lower_thres_intensity))
    
    # Colorize if needed
    if args.transparent_below > 0:
        mask = (chunk < args.transparent_below/255)
    if args.colorize is not None:
        if args.transparent_below > 0:
            chunk = np.stack([chunk * r, chunk * g, chunk * b, chunk * 255], axis=-1)
            chunk[mask, 3] = 0 ## make it transparent
        else:
            chunk = np.stack([chunk * r, chunk * g, chunk * b], axis=-1)
    else:
        if args.transparent_below > 0:
            chunk = np.stack([chunk * 255, chunk * 255, chunk * 255, chunk * 255], axis=-1)
            chunk[mask, 3] = 0 ## make it transparent
        else:
            chunk = chunk * 255.0
    
    return chunk.astype(np.uint8)

def compute_quantiles_from_histogram(histogram, total_count, quantile):
    """
    Given a histogram (a mapping from integer values to counts) and the total count of values,
    compute the integer corresponding to each quantile in the list `quantiles`.
    
    Args:
        histogram: Counter mapping integer values to counts.
        total_count: Total number of elements processed.
        quantiles: A list of quantile levels (floats between 0 and 1).
        
    Returns:
        A dictionary mapping each quantile level to the corresponding integer value.
    """
    # Sort the keys (the integer values) in increasing order
    sorted_keys = sorted(histogram.keys())
    # Create a cumulative count array
    cum_counts = []
    running_total = 0
    for key in sorted_keys:
        running_total += histogram[key]
        cum_counts.append(running_total)
    
    # For each desired quantile, find the corresponding integer value
    quantile_values = {}
    target = quantile * total_count
        
    # Find the first index where the cumulative count is >= target.
    # np.searchsorted does exactly that.
    idx = np.searchsorted(cum_counts, target, side='left')
    # The corresponding integer value is:
    return sorted_keys[idx]

def image_ome2png(_args):
    args = parse_arguments(_args)
    logger = create_custom_logger(__name__, args.out_prefix + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    out_dir=os.path.dirname(args.out_prefix)
    os.makedirs(out_dir, exist_ok=True)
    
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

    # [CSV processing and metadata extraction remain the same]
    with tifffile.TiffFile(args.tif, _multifile=False) as tif:
        logger.info(f"Loaded OME-TIFF file {args.tif}")
        
        n_pages = len(tif.pages)
        n_series = len(tif.series)
        n_levels = len(tif.series[0].levels)
        
        if n_pages > 1 and args.page is None:
            logger.error(f"Multiple pages detected (n={n_pages}). Please specify the page number to extract the image from")
            sys.exit(1)

        # Get page and validate
        args.page = 0 if args.page is None else args.page
        page = tif.series[args.series].levels[args.level].pages[args.page]
        
        #assert page.is_tiled, "Only tiled TIFF files are supported"
        if not args.high_memory and not page.is_tiled:
            logger.error("When the TIFF file is not tiled, please use the --high-memory flag")
            sys.exit(1)
        
        if len(page.shape) == 3:
            if page.shape[2] != 3:
                logger.error("The colored image is not in RGB format")
                sys.exit(1)
            is_mono = False
            if args.colorize:
                logger.error("Cannot colorize already RGB-colored image")
                sys.exit(1)
        elif args.colorize is not None:
            is_mono = False
        else:
            is_mono = True

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
        
        ## iterate over the segment and extract the image by tile
        if args.high_memory:
            (chunk_height, chunk_width) = page.shape[:2]
            n_chunks = 1
        else:
            (chunk_height, chunk_width) = page.chunks
            n_chunks = ((page.imagelength + chunk_height - 1) // chunk_height) * ((page.imagewidth + chunk_width - 1) // chunk_width)
        logger.info(f"Processing in chunks of {chunk_height}x{chunk_width} pixels")

        #pixel_bytes = 4 if args.colorize else 1  # Account for RGB vs grayscale
            
        # Prepare output file
        output_shape = (*page.shape[:2], 4) if args.transparent_below > 0 else ((*page.shape[:2], 3) if args.colorize else page.shape[:2])
        output_dtype = np.uint8
        
        if args.high_memory:
            image_array_highmem = page.asarray()
        
        ## if needed, get the overall distribution of pixel values for quantile-based thresholding
        if args.upper_thres_quantile is not None or args.lower_thres_quantile is not None:
            logger.info("Calculating pixel value distribution...")
            histogram = Counter()
            total_count = 0
            if args.high_memory:
                segments = [image_array_highmem]
            else:
                segments = page.segments()

            for i, segment in enumerate(segments):
                if i % 100 == 0:
                    logger.info(f"Processing segment {i+1}/{n_chunks}...")
                data = np.squeeze(segment[0])
                flat = data.ravel()
                histogram.update(flat.tolist())
                total_count += flat.size
            
            if args.upper_thres_quantile is not None:
                #args.upper_thres_intensity = np.quantile(pixel_values, args.upper_thres_quantile)
                args.upper_thres_intensity = compute_quantiles_from_histogram(histogram, total_count, args.upper_thres_quantile)
                logger.info(f"Upper threshold for quantile {args.upper_thres_quantile} is set to {args.upper_thres_intensity}")
                
            if args.lower_thres_quantile is not None:
                #args.lower_thres_intensity = np.quantile(pixel_values, args.lower_thres_quantile)
                args.lower_thres_intensity = compute_quantiles_from_histogram(histogram, total_count, args.lower_thres_quantile)
                logger.info(f"Lower threshold for quantile {args.lower_thres_quantile} is set to {args.lower_thres_intensity}")
            del histogram
            gc.collect()
    
        # Create memory-mapped output file
        output = np.memmap(f"{args.out_prefix}_output.npy", dtype=output_dtype, mode='w+', shape=output_shape)
    
        # Process color information
        r, g, b = None, None, None
        if args.colorize:
            r, g, b = hex_to_rgb(args.colorize)


        if args.high_memory:
            segments = [image_array_highmem]
        else:
            segments = page.segments()
            
        for i, segment in enumerate(segments):
            if i % 100 == 0:
                current_memory = get_memory_usage()
                logger.info(f"Processing segment {i+1}/{n_chunks}, using {current_memory:.2f} GB memory...")

            # [Segment processing code remains the same]
            (data, offset, bytecount) = segment
            offset_y, offset_x = offset[-3], offset[-2]
            
            ## determine height and length
            if offset_y + chunk_height > page.imagelength:
                height = page.imagelength - offset_y
            else:
                height = chunk_height
                
            if offset_x + chunk_width > page.imagewidth:
                width = page.imagewidth - offset_x
            else:
                width = chunk_width

            data = np.squeeze(data)
            processed = process_image_chunk(data, args, r, g, b)
            
            #print(f"Chunk {i}: shape: {data.shape}, {processed.shape}, {offset_x}, {offset_y}, {width}, {height}, {page.imagewidth}, {page.imagelength}")
                
            # Write to output
            if len(output_shape) > 2:
                output[offset_y:(offset_y+height), offset_x:(offset_x+width), :] = processed[0:height, 0:width, :]
            else:
                output[offset_y:(offset_y+height), offset_x:(offset_x+width)] = processed[0:height, 0:width]
            
            del data, processed, segment
            #gc.collect()
            #output.flush()            
                
        # Save final image
        output.flush()
        logger.info("Saving final image...")
        Image.fromarray(output).save(f"{args.out_prefix}.png")
        
        # save the bounds
        with open(f"{args.out_prefix}.bounds.csv", 'w') as f:
            f.write(f"{ul[0]},{ul[1]},{lr[0]},{lr[1]}\n")

        # save the color mode 
        if args.write_color_mode:
            color_mode = "rgb"
            if args.transparent_below > 0:
                color_mode = "rgba"
            elif is_mono:
                color_mode = "mono"
            with open(f"{args.out_prefix}.color.csv", 'w') as f:
                f.write(f"{color_mode}\n")

        # Clean up temporary files
        os.remove(f"{args.out_prefix}_output.npy")
        if os.path.exists(f"{args.out_prefix}_transformed.npy"):
            os.remove(f"{args.out_prefix}_transformed.npy")

    logger.info("Analysis Completed")

if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], script_name)
    func(sys.argv[1:])
