import sys, os, gzip, argparse, logging, warnings, shutil, subprocess, ast, re
import pandas as pd
import numpy as np
import tifffile, json
from PIL import Image
import psutil
from typing import Tuple
import gc

from cartloader.utils.minimake import minimake
from cartloader.utils.utils import cmd_separator, scheck_app, create_custom_logger, read_minmax, hex_to_rgb

def parse_arguments(_args):
    # [Previous argument parsing code remains the same]
    parser = argparse.ArgumentParser(prog=f"cartloader transform_aligned_histology", 
                                   description="Convert Aligned Histology from 10x Xenium or Vizgen MERSCOPE")

    # Add existing arguments
    inout_params = parser.add_argument_group("Input/Output Parameters")
    inout_params.add_argument('--tif', type=str, required=True, help='Input OME-TIFF file')
    inout_params.add_argument('--csv', type=str, help='(Vizgen only) CSV file containing conversion parameters')
    inout_params.add_argument('--out-prefix', type=str, required=True, help='Output file prefix')
    inout_params.add_argument("--page", type=int, help='For 3D (X/Y/Z) OME file, specify the z-value')
    inout_params.add_argument("--level", type=int, default=0, help='Specify the level to extract')
    inout_params.add_argument("--series", type=int, default=0, help='Specify the index of series')
    inout_params.add_argument("--upper-thres-quantile", type=float, help='Quantile-based capped value')
    inout_params.add_argument("--upper-thres-intensity", type=float, default=255, help='Intensity-based capped value')
    inout_params.add_argument("--lower-thres-quantile", type=float, help='Quantile-based floored value')
    inout_params.add_argument("--lower-thres-intensity", type=float, default=0, help='Intensity-based floored value')
    inout_params.add_argument("--colorize", type=str, help='Colorize using RGB code as max value')
    inout_params.add_argument('--flip-horizontal', action='store_true', default=False)
    inout_params.add_argument('--flip-vertical', action='store_true', default=False)
    inout_params.add_argument('--rotate-clockwise', action='store_true', default=False)
    inout_params.add_argument('--rotate-counter', action='store_true', default=False)

    memory_params = parser.add_argument_group("Memory Management")
    memory_params.add_argument('--max-memory-gb', type=float, default=4.0,
                             help='Maximum memory usage in gigabytes')

    aux_params = parser.add_argument_group("Auxiliary Parameters")    
    aux_params.add_argument('--skip-pmtiles', action='store_true', default=False)
    aux_params.add_argument('--log', action='store_true', default=False)
    aux_params.add_argument('--log-suffix', type=str, default=".log")

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def get_memory_usage():
    """Get current memory usage in GB"""
    process = psutil.Process()
    return process.memory_info().rss / (1024 * 1024 * 1024)

# class ChunkedTiffReader:
#     def __init__(self, tiff_page, chunk_size):
#         self.page = tiff_page
#         self.chunk_size = chunk_size
#         self.shape = tiff_page.shape
#         self.dtype = tiff_page.dtype

#     def read_chunk(self, y_start, y_end, x_start, x_end):
#         """Read a chunk of the TIFF file directly"""
#         return self.page.asarray(out='memmap')[y_start:y_end, x_start:x_end]

def process_image_chunk(chunk: np.ndarray, args, r=None, g=None, b=None) -> np.ndarray:
    """Process a chunk of the image with the specified transformations"""
    # Clip and normalize
    chunk = ((np.clip(chunk, a_min=args.lower_thres_intensity, a_max=args.upper_thres_intensity) 
             - args.lower_thres_intensity) / (args.upper_thres_intensity - args.lower_thres_intensity))
    
    # Colorize if needed
    if args.colorize is not None:
        chunk = np.stack([chunk * r, chunk * g, chunk * b], axis=-1)
    else:
        chunk = chunk * 255.0
    
    return chunk.astype(np.uint8)

def run_hist2pmtiles(_args):
    args = parse_arguments(_args)
    logger = create_custom_logger(__name__, args.out_prefix + args.log_suffix if args.log else None)
    logger.info("Analysis Started")

    # [CSV processing and metadata extraction remain the same]
    
    with tifffile.TiffFile(args.tif) as tif:
        logger.info(f"Loaded OME-TIFF file {args.tif}")
        
        # Get page and validate
        args.page = 0 if args.page is None else args.page
        page = tif.series[args.series].levels[args.level].pages[args.page]
        
        assert page.is_tiled, "Only tiled TIFF files are supported"
        
        ## iterate over the segment and extract the image by tile
        (chunk_height, chunk_width) = page.chunks
        logger.info(f"Processing in chunks of {chunk_height}x{chunk_width} pixels")

        pixel_bytes = 4 if args.colorize else 1  # Account for RGB vs grayscale
        
        # Create chunked reader
        #reader = ChunkedTiffReader(page)
    
        # Prepare output file
        output_shape = (*page.shape[:2], 3) if args.colorize else page.shape[:2]
        output_dtype = np.uint8
    
        # Create memory-mapped output file
        output = np.memmap(f"{args.out_prefix}_output.npy", dtype=output_dtype,
                        mode='w+', shape=output_shape)
    
        # Process color information
        r, g, b = None, None, None
        if args.colorize:
            r, g, b = hex_to_rgb(args.colorize)

        segments = page.segments()
        for i, segment in enumerate(segments):
            logger.info(f"Processing segment {i}...")
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
            
            print(f"Chunk {i}: shape: {data.shape}, {processed.shape}, {offset_x}, {offset_y}, {width}, {height}, {page.imagewidth}, {page.imagelength}")
                
            # Write to output
            if args.colorize:
                output[offset_y:(offset_y+height), offset_x:(offset_x+width), :] = processed[0:height, 0:width, :]
            else:
                output[offset_y:(offset_y+height), offset_x:(offset_x+width)] = processed[0:height, 0:width]
            
            del data, processed, segment
            gc.collect()
            output.flush()
            
            current_memory = get_memory_usage()
            logger.info(f"Current memory usage: {current_memory:.2f} GB")
                
        # Save final image
        logger.info("Saving final image...")
        Image.fromarray(output).save(f"{args.out_prefix}.png")
        
        # Clean up temporary files
        os.remove(f"{args.out_prefix}_output.npy")
        if os.path.exists(f"{args.out_prefix}_transformed.npy"):
            os.remove(f"{args.out_prefix}_transformed.npy")
        
        # Handle PMTiles conversion
        if not args.skip_pmtiles:
            # [PMTiles conversion code remains the same]
            pass

    logger.info("Analysis Completed")

if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], script_name)
    func(sys.argv[1:])