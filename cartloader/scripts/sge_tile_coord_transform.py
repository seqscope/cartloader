import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess
import numpy as np

from cartloader.utils.utils import log_dataframe, create_custom_logger

def parse_minmax(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Minmax file not found: {file_path}")
    with open(file_path, "r") as f:
        return {k: float(v) for k, v in (line.strip().split("\t") for line in f)}

# This is separated from combine_sge_by_layout.py to be reusable for hist_stitch process.
def sge_tile_coord_transform(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Calculate the offsets and global minmax for each tile SGE"    
    )
    parser.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of input information in a specific format: <minmax_path>,<row>,<col>.")
    parser.add_argument('--output', type=str, help='Output subset file name for the coordinate minmax and offset TSV file .') #"coordinate_minmax_per_tile.tsv"
    parser.add_argument("--units-per-um", type=float, default=1.0, help="If --convert-to-um, define the scaling factor for unit conversion from coordinate to um(default: 1.0)")
    parser.add_argument('--log', action='store_true', default=False, help='Write log to file')
    args = parser.parse_args(_args)

    # outdir
    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True) 

    # log
    logger = create_custom_logger(__name__, os.path.join(out_dir, f"{args.output}.log") if args.log else None)
    logger.info(f"Command: {' '.join(sys.argv)}")

    # 2. input
    # * input paths
    assert len(args.in_tiles) > 0, "No input tiles provided."
    assert all(len(in_tile.split(",")) == 3 for in_tile in args.in_tiles), "Each input tile should have 3 elements: <minmax_path>,<row>,<col>."

    df = pd.DataFrame(args.in_tiles, columns=["input"])
    df = df["input"].str.split(",", expand=True)
    df.columns = ["minmax_path", "row", "col"]
    df["row"] = df["row"].astype(int)
    df["col"] = df["col"].astype(int)

    # * input min max 
    minmax_data = df["minmax_path"].apply(parse_minmax)
    minmax_df = pd.DataFrame(minmax_data.tolist())
    assert {"xmin", "xmax", "ymin", "ymax"}.issubset(minmax_df.columns), "Missing required minmax columns: xmin, xmax, ymin, ymax in the minmax file."

    # * concat
    df = pd.concat([df, minmax_df], axis=1)
    log_dataframe(df, msg="  - Input SGEs with local coords:")

    # 3. Calculate subset (offset and minmax for each tile)
    logger.info(f"  "+"-"*20)
    logger.info(f"  - Calculating offsets and new minmax for each tile ...")
    # (1) The tile should be put at the center of the fixed space.

    # * All x and y across rows and cols (in unit)
    df["tile_x"] = df["xmax"] - df["xmin"]  # tile width (x)
    df["tile_y"] = df["ymax"] - df["ymin"]  # tile height (y)

    # * Fixed tile size (using the largest x and y) per row and per col (in unit)
    row2y = df.groupby("row")["tile_y"].max() 
    col2x = df.groupby("col")["tile_x"].max()
    logger.info(f"      * Fixed height per row (in coordinate unit): {row2y.to_dict()}")
    logger.info(f"      * Fixed width per col (in coordinate unit): {col2x.to_dict()}")

    # * Offset for each col and per row (in unit)
    row2y_offsets = row2y.cumsum() - row2y / 2
    col2x_offsets = col2x.cumsum() - col2x / 2

    # * drop the tile x and y
    df.drop(["tile_x", "tile_y"], axis=1, inplace=True)

    # * Offset for each tile (in unit)
    df["x_offset"] = df["col"].map(col2x_offsets) - (df["xmin"] + df["xmax"]) / 2
    df["y_offset"] = df["row"].map(row2y_offsets) - (df["ymin"] + df["ymax"]) / 2

    df_subset = df[["row", "col", "xmin", "xmax", "ymin", "ymax", "x_offset", "y_offset"]].copy()

    df_subset["units_per_um"]= args.units_per_um

    # use a scaled version (only for record)
    df_subset["global_xmin_um"] = (df_subset["xmin"] + df_subset["x_offset"])/ args.units_per_um
    df_subset["global_xmax_um"] = (df_subset["xmax"] + df_subset["x_offset"])/ args.units_per_um
    df_subset["global_ymin_um"] = (df_subset["ymin"] + df_subset["y_offset"])/ args.units_per_um
    df_subset["global_ymax_um"] = (df_subset["ymax"] + df_subset["y_offset"])/ args.units_per_um

    df_subset.rename(columns={ "xmin": "local_xmin_unit",
                               "xmax": "local_xmax_unit",
                                "ymin": "local_ymin_unit",
                                "ymax": "local_ymax_unit",
                                "x_offset": "x_offset_unit",
                                "y_offset": "y_offset_unit"}, inplace=True)

    df_subset.to_csv(args.output, sep="\t", index=False)
    logger.info(f"  - Offsets and new minmax per tile written to {args.output}")

if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
