import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess
import numpy as np
from collections import defaultdict

from cartloader.utils.utils import log_dataframe, create_custom_logger

def parse_minmax(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Minmax file not found: {file_path}")
    with open(file_path, "r") as f:
        return {k: float(v) for k, v in (line.strip().split("\t") for line in f)}

# This is separated from combine_sge_by_layout.py to be reusable for hist_stitch process.
def feature_overlapping(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="""
        Identify features shared across all input tiles and generate summary outputs.
        - The input tiles (--in-tiles) are always required and used to determine which features are shared.
        - Optionally, a separate feature file (--in-feature) can be supplied to extract summary count data instead of pulling this from the tiles directly.
        - If --min-ct-per-ftr-tile > 0, features are filtered based on a minimum count threshold across all tiles.
        """
        
    )
    parser.add_argument("--in-tiles", type=str, nargs='*',  required=True, default=[], help="List of input information in a specific format: <feature_path>,<row>,<col>.")
    parser.add_argument('--in-feature', type=str, default=None, help='(Optional) Optional. File with feature names and counts. Used to populate --out-overlap if provided. (default: None)')
    parser.add_argument('--out-overlap', type=str, required=True, help='Output file for only the features that are present in all tiles. The output columns includes <feature_name>, <sum of count across all tiles> per count.')
    parser.add_argument('--out-dist', type=str, help='(Optional) Output distribution file for all features. The output columns includes <feature_name>, <number of tiles with the feature>, and <count> per tile per count.')
    parser.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    parser.add_argument("--colnames-count", type=str, nargs='*', help="Columns (default: count).", default=['count'])
    parser.add_argument('--min-ct-per-ftr-tile', type=int, default=0, help='Minimum count to keep a shared feature (default: 0).')
    parser.add_argument('--colname-key-count', type=str, default=None, help='Column to use for count filtering. Defaults to first --colnames-count if only one')
    parser.add_argument('--log', action='store_true', default=False, help='Write log to file')
    args = parser.parse_args(_args)

    # outdir
    out_dir = os.path.dirname(args.out_overlap)
    os.makedirs(out_dir, exist_ok=True) 

    # log
    logger = create_custom_logger(__name__, os.path.join(out_dir, f"{args.out_dist}.log") if args.log else None)
    logger.info(f"Command: {' '.join(sys.argv)}")

    # 2. input
    # * input paths
    num_of_tiles = len(args.in_tiles)
    assert num_of_tiles > 0, "No input tiles provided."
    assert all(len(in_tile.split(",")) == 3 for in_tile in args.in_tiles), "Each input tile should have 3 elements: <minmax_path>,<row>,<col>."

    df = pd.DataFrame(args.in_tiles, columns=["input"])
    df = df["input"].str.split(",", expand=True)
    df.columns = ["feature_path", "row", "col"]
    df["row"] = df["row"].astype(int)
    df["col"] = df["col"].astype(int)

    # Start with an empty DataFrame
    ftr_data = None

    for idx, row in df.iterrows():
        feature_path = row["feature_path"]
        row_id = row["row"]
        col_id = row["col"]

        logger.info(f"Reading {feature_path} for row {row_id} col {col_id}")
        data = pd.read_csv(feature_path, sep='\\t')

        # Select and rename relevant columns
        selected = data[[args.colname_feature_name ] + args.colnames_count].copy()
        selected.set_index(args.colname_feature_name, inplace=True)
        selected.rename(columns={c: f"{row_id}_{col_id}_{c}" for c in args.colnames_count}, inplace=True)

        # Merge with running ftr_data
        if ftr_data is None:
            ftr_data = selected
        else:
            ftr_data = ftr_data.merge(selected, how='outer', left_index=True, right_index=True)

    # feature distribution
    ftr_data.insert(0, "N", ftr_data.notna().sum(axis=1))
    ftr_data.insert(0, args.colname_feature_name, ftr_data.index)
    ftr_data.reset_index(drop=True, inplace=True)

    if args.out_dist:
        ftr_data.to_csv(args.out_dist, sep="\t", index=False, na_rep="NA", compression="gzip")
        logger.info(f"Feature distribution saved to {args.out_dist}")

    # filter1: keep features existing in all tiles
    logger.info(f"* Number of input features: {ftr_data.shape[0]}")
    overlap_ftr_data = ftr_data[ftr_data["N"] == num_of_tiles]
    logger.info(f"* Number of features shared across all tiles: {overlap_ftr_data.shape[0]}")

    # filter2: keep features min count across all tiles > threshold
    if args.min_ct_per_ftr_tile > 0:
        if args.colname_key_count is None and len(args.colnames_count) == 1:
            args.colname_key_count = args.colnames_count[0]
        assert args.colname_key_count is not None, "Please specify --colname-key-count or use --colnames-count with only one column when filtering by min count."
        
        key_counts = [f"{row_id}_{col_id}_{args.colname_key_count}" for row_id, col_id in zip(df["row"], df["col"])]
        overlap_ftr_data["min_count"] = overlap_ftr_data[key_counts].min(axis=1)
        
        overlap_ftr_data = overlap_ftr_data[overlap_ftr_data["min_count"] >= args.min_ct_per_ftr_tile]
        logger.info(f"* Number of shared features with a minimum count of {args.min_ct_per_ftr_tile} across all tiles: {overlap_ftr_data.shape[0]}")

    # extract feature names
    overlap_ftrs = overlap_ftr_data[args.colname_feature_name].tolist()

    # extract feature and counts
    if args.in_feature:
        df_ftr = pd.read_csv(args.in_feature, sep="\t")
        logger.info(f"Reading a separate feature file {args.in_feature} with {df_ftr.shape[0]} features.")
        df_ftr = df_ftr[df_ftr[args.colname_feature_name].isin(overlap_ftrs)]
    else:
        logger.info(f"Using the features from the input tiles.")
        for col_count in args.colnames_count:
            overlap_ftr_data[col_count] = overlap_ftr_data[[f"{row_id}_{col_id}_{col_count}" for row_id, col_id in zip(df["row"], df["col"])]].sum(axis=1)
        df_ftr = overlap_ftr_data

    logger.info(f"Returning {df_ftr.shape[0]} features with counts.")
    df_ftr = df_ftr[[args.colname_feature_name] + args.colnames_count]
    for col_count in args.colnames_count:
        df_ftr[col_count] = df_ftr[col_count].astype(int)
    df_ftr.to_csv(args.out_overlap, sep="\t", index=False, compression="gzip")
    logger.info(f"Feature overlap saved to {args.out_overlap}")

if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
