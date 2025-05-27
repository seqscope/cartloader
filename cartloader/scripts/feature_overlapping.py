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
        description=""
    )
    parser.add_argument("--in-dist", type=str, required=True, help="Input feature distribution file")
    parser.add_argument('--in-feature', type=str, default=None, required=True, help='Input feature file')
    parser.add_argument('--output', type=str, required=True, help='Output file for only the features that are present in all tiles. The output columns includes <feature_name>, <sum of count across all tiles> per count.')
    parser.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    parser.add_argument("--colname-count", type=str,  default="count", help="Column name for the count (default: count).")
    parser.add_argument("--colname-n", type=str, default='N', help="Column name for the number of tiles (default: N).")
    parser.add_argument('--min-ct-per-ftr-tile', type=int, default=0, help='Minimum count to keep a shared feature (default: 0).')
    parser.add_argument('--log', action='store_true', default=False, help='Write log to file')
    args = parser.parse_args(_args)

    # outdir
    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True) 

    # log
    logger = create_custom_logger(__name__, os.path.join(out_dir, f"{args.output}.log") if args.log else None)
    logger.info(f"Command: {' '.join(sys.argv)}")

    # 2. locate features
    df = pd.read_csv(args.in_dist, sep="\t", header=0)
    
    tilecount_colnames = [c for c in df.columns.tolist() if c not in [args.colname_feature_name, args.colname_n]]
    tilecount_colnames = [c for c in tilecount_colnames if re.match(r"^\d+_\d+_" + args.colname_count, str(c))]
    assert len(tilecount_colnames) > 0, f"Cannot find any tile count columns in {args.in_dist}."
    logger.info(f"Tile count columns: {tilecount_colnames}")

    n_tiles = len(tilecount_colnames)
    logger.info(f"Number of tiles: {n_tiles}")

    logger.info(f"Feature:")
    logger.info(f"* Number of input features: {df.shape[0]}")
    print(args.colname_n)
    print(n_tiles)
    df = df[df[args.colname_n] == n_tiles]
    logger.info(f"* Number of features shared across all tiles: {df.shape[0]}")

    # filter2: keep features min count across all tiles > threshold
    df = df[[args.colname_feature_name] + tilecount_colnames]
    df["min_count"] = df[tilecount_colnames].min(axis=1)
    df = df[df["min_count"] >= args.min_ct_per_ftr_tile]
    logger.info(f"* Number of shared features with a minimum count of {args.min_ct_per_ftr_tile} across all tiles: {df.shape[0]}")
    overlap_ftrs= df[args.colname_feature_name].tolist()

    # 3. read the feature file
    df_ftr = pd.read_csv(args.in_feature, sep="\t")
    logger.info(f"Input feature: {args.in_feature}")
    logger.info(f" * Input N feature: {df_ftr.shape[0]} features.")
    df_ftr = df_ftr[df_ftr[args.colname_feature_name].isin(overlap_ftrs)]
   
    logger.info(f" * After filtering: {df_ftr.shape[0]} features.")
    df_ftr[args.colname_count] = df_ftr[args.colname_count].astype(int)
    df_ftr.to_csv(args.output, sep="\t", index=False, compression="gzip")
    logger.info(f"Feature overlap saved to {args.output}")

if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
