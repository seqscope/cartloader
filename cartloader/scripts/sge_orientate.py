import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
from collections import defaultdict
import math
from cartloader.utils.image_helper import update_orient, orient2axisorder

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""Orientate SGE""")
    parser.add_argument("--in-transcript", type=str, help="Input file.", required=True)
    parser.add_argument("--out-transcript", type=str, help="Output file.", required=True)
    parser.add_argument("--in-feature", type=str, help="Input feature file.", required=True)
    parser.add_argument("--out-feature", type=str, help="Output feature file.", required=True)
    parser.add_argument("--in-minmax", type=str, help="Input minmax file.", required=True)
    parser.add_argument("--out-minmax", type=str, help="Output minmax file.", required=True)
    parser.add_argument("--chunk-size", type=int, default=10**6, help="Chunk size for processing (default: 10^6).")
    parser.add_argument('--rotate', type=str, default=None, choices=["90", "180", "270"],  help='Rotate by 90, 180, or 270 degrees clockwise. Rotate precedes flip.')
    parser.add_argument('--flip-vertical', action='store_true', default=False, help='Flip vertically (flipped along Y-axis). Rotate precedes flip.')
    parser.add_argument('--flip-horizontal', action='store_true', default=False, help='Flip horizontally (flipped along X-axis). Rotate precedes flip.')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)

def orient_minmax(minmax, axis_order):
    
    order_x, order_y = axis_order.strip().split(",")

    def get_bounds(sign, minval, maxval):
        if sign == "-":
            return (
                -maxval + minval + maxval,
                -minval + minval + maxval,
            )
        else:
            return (minval, maxval)

    bounds = {
        "1": get_bounds("", minmax["xmin"], minmax["xmax"]),
        "-1": get_bounds("-", minmax["xmin"], minmax["xmax"]),
        "2": get_bounds("", minmax["ymin"], minmax["ymax"]),
        "-2": get_bounds("-", minmax["ymin"], minmax["ymax"]),
    }

    return {
        "xmin": min(bounds[order_x]),
        "xmax": max(bounds[order_x]),
        "ymin": min(bounds[order_y]),
        "ymax": max(bounds[order_y]),
    }

def orient_transcript_chunk(df_chunk, minmax, axis_order):
    order_x, order_y = axis_order.strip().split(",")

    # Helper to fetch and flip coordinate
    def get_transformed(coord, sign, minval, maxval):
        return -coord + maxval + minval if sign == "-" else coord

    # Create new coordinates
    coord_map = {
        "1": get_transformed(df_chunk["X"], "", minmax["xmin"], minmax["xmax"]),
        "-1": get_transformed(df_chunk["X"], "-", minmax["xmin"], minmax["xmax"]),
        "2": get_transformed(df_chunk["Y"], "", minmax["ymin"], minmax["ymax"]),
        "-2": get_transformed(df_chunk["Y"], "-", minmax["ymin"], minmax["ymax"]),
    }

    df_chunk = df_chunk.copy()
    df_chunk["X"] = coord_map[order_x]
    df_chunk["Y"] = coord_map[order_y]
    return df_chunk

def sge_orientate(_args):
    args = parse_arguments(_args)

    assert os.path.exists(args.in_transcript), f"Input file {args.in_transcript} does not exist."
    assert os.path.exists(args.in_minmax), f"Input minmax file {args.in_minmax} does not exist."

    rotate, flip_vertical, flip_horizontal = update_orient (args.rotate, args.flip_vertical, args.flip_horizontal, f"sge:{args.in_transcript}")
    axis_order = orient2axisorder.get((rotate, flip_vertical, flip_horizontal))
    
    # orientate minmax
    minmax_data= pd.read_csv(args.in_minmax, sep="\t", header=None)
    minmax_data[1] = minmax_data[1].astype(float)
    minmax = dict(zip(minmax_data[0], minmax_data[1]))

    new_minmax = orient_minmax(minmax, axis_order)

    with open(args.out_minmax, "w") as f:
        for key, value in new_minmax.items():
            f.write(f"{key}\t{value}\n")
    print(f"Orientated minmax file saved to {args.out_minmax}")

    # copy the feature file
    shutil.copy(args.in_feature, args.out_feature)
    print(f"Feature file copied to {args.out_feature}")

    # orientate transcript
    first_chunk = True
    for chunk in pd.read_csv(args.in_transcript, sep="\t", chunksize=args.chunk_size):
        transformed = orient_transcript_chunk(chunk, minmax, axis_order)
        transformed.to_csv(
            args.out_transcript,
            sep="\t",
            index=False,
            header=first_chunk,
            mode='wt' if first_chunk else 'at',
            compression="gzip"
        )
        first_chunk = False
    
    print(f"Orientated transcript file saved to {args.out_transcript}")

if __name__ == "__main__":
    # get the cartloader path
    global cartloader_repo
    cartloader_repo=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])