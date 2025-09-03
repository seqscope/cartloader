import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess
import rasterio # for extracting bounds

def image_extract_ullr(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="""
        Extracts the ULLR (Upper Left, Lower Right) coordinates from a TIFF file.
        Outputs: A text file with one line of Upper Left X, Upper Left Y, Lower Right X, Lower Right Y.
        Options: Applying X/Y offsets or cropping using specified X Y min max bounds. 
        """
    )
    parser.add_argument("--input", type=str, help="Input tif file.", required=True)
    parser.add_argument("--output", type=str, help="Output file.", required=True)
    parser.add_argument("--x-offset", type=float, default=0, help="Value to add to the X coordinates (default: 0).")
    parser.add_argument("--y-offset", type=float, default=0, help="Value to add to the Y coordinates (default: 0).")
    parser.add_argument("--crop-by-minmax", action="store_true", help="Enable cropping using --minmax bounds.")
    parser.add_argument("--minmax", type=str, help="Comma-separated list defining the cropping bounds as xmin,xmax,ymin,ymax.")
    args= parser.parse_args(_args)
    
    with rasterio.open(args.input) as src:
        bounds=src.bounds
        transform = src.transform # a 2D affine transformation matrix. 
        
    if not args.crop_by_minmax:
        # Use native image bounds
        ulx = bounds.left + args.x_offset
        uly = bounds.top + args.y_offset
        lrx = bounds.right + args.x_offset
        lry = bounds.bottom + args.y_offset
    else:
        assert args.minmax, "--minmax is required when --crop-by-minmax is set"
        xmin, xmax, ymin, ymax = map(float, args.minmax.split(","))

        # Determine orientation using transform
        pixel_width = transform.a
        pixel_height = transform.e

        # ULLR assignment based on direction of axes
        if pixel_width > 0 and pixel_height < 0:
            ulx, uly = xmin, ymax
            lrx, lry = xmax, ymin
        elif pixel_width < 0 and pixel_height < 0:
            ulx, uly = xmax, ymax
            lrx, lry = xmin, ymin
        elif pixel_width > 0 and pixel_height > 0:
            ulx, uly = xmin, ymin
            lrx, lry = xmax, ymax
        elif pixel_width < 0 and pixel_height > 0:
            ulx, uly = xmax, ymin
            lrx, lry = xmin, ymax
        else:
            raise ValueError("Unsupported image orientation (affine): only axis-aligned (no rotation/shear) with non-zero scales are supported; reproject to north-up first")

        ulx += args.x_offset
        uly += args.y_offset
        lrx += args.x_offset
        lry += args.y_offset

    with open(args.output, "w") as f:
        f.write(f"{ulx} {uly} {lrx} {lry}\n")

    print(f"Output written to {args.output}")
if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
