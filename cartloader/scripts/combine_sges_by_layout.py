import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess
import numpy as np

from cartloader.utils.utils import log_dataframe, create_custom_logger

global script_name
script_name = os.path.splitext(os.path.basename(__file__))[0]

def parse_minmax(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Minmax file not found: {file_path}")
    with open(file_path, "r") as f:
        return {k: float(v) for k, v in (line.strip().split("\t") for line in f)}

def combine_sges_by_layout(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Combine SGE data by layout. Note the transcript/feature/minmax files should be in the same resolution.",
    )
    parser.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of input information in a specific format: <transcript_path>,<feature_path>,<minmax_path>,<row>,<col>.")
    parser.add_argument("--out-dir", type=str, help="Output directory.")
    parser.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='Output file name for the compressed transcript-indexed SGE file in TSV format (default: transcripts.unsorted.tsv.gz).')
    parser.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='Output minmax file name for the coordinate minmax TSV file (default: coordinate_minmax.tsv).')
    parser.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='Output feature file name for the compressed UMI count per gene TSV file (default: feature.clean.tsv.gz).')
    parser.add_argument('--out-subset', type=str, default="subset_minmax.tsv", help='Output subset file name for the coordinate minmax and offset TSV file (default: subset_minmax.tsv).')
    parser.add_argument("--colnames-count", type=str, nargs='*', help="Columns to sum (default: count).", default=['count'])
    parser.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    parser.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    parser.add_argument('--colname-x', type=str, default="X", help='X column name (default: X)')
    parser.add_argument('--colname-y', type=str, default="Y", help='Y column name (default: Y)')
    parser.add_argument("--convert-to-um", action="store_true", help="Convert output coordinates to micrometers(default: False)")
    parser.add_argument("--units-per-um", type=float, default=1.0, help="If --convert-to-um, define the scaling factor for unit conversion from coordinate to um(default: 1.0)")
    # parser.add_argument('--minmax-in-um', action='store_true', help='Input minmax is in um while input transcript is based on unit.') # This was originally designed for the minmax in SeqScope v2 SGEs. Now we only take the transcripts and generate minmax from the sge-adds-on. So, disabling this.
    parser.add_argument("--precision", type=int, default=2, help="Precision for the output minmax and transcript files (default: 2)")
    parser.add_argument('--debug', action='store_true', help='Test mode.')
    args = parser.parse_args(_args)

    # outdir
    if os.path.isfile(args.out_dir):
        raise NotADirectoryError(f"Output path exists as a file, not a directory: {args.out_dir}")

    os.makedirs(args.out_dir, exist_ok=True) 

    # log
    log_f = os.path.join(args.out_dir, f"sge_stitch.log")
    logger = create_custom_logger(__name__, log_f)
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"  - Output directory: {args.out_dir}")

    # 1. output cols
    out_cols=[args.colname_x, args.colname_y]
    out_cols.extend([args.colname_feature_name]) if args.colname_feature_id is None else out_cols.extend([args.colname_feature_name, args.colname_feature_id])
    out_cols.extend(args.colnames_count)
    logger.info(f"  - Output columns: {out_cols}")

    # 2. input
    # * input paths
    assert len(args.in_tiles) > 0, "No input tiles provided."
    assert all(len(in_tile.split(",")) == 5 for in_tile in args.in_tiles), "Each input tile should have 5 elements: <transcript_path>,<feature_path>,<minmax_path>,<row>,<col>."

    df = pd.DataFrame(args.in_tiles, columns=["input"])
    df = df["input"].str.split(",", expand=True)
    df.columns = ["transcript_path", "feature_path", "minmax_path", "row", "col"]
    df["row"] = df["row"].astype(int)
    df["col"] = df["col"].astype(int)

    # * input min max 
    minmax_data = df["minmax_path"].apply(parse_minmax)
    minmax_df = pd.DataFrame(minmax_data.tolist())
    assert {"xmin", "xmax", "ymin", "ymax"}.issubset(minmax_df.columns), "Missing required minmax columns: xmin, xmax, ymin, ymax in the minmax file."

    # * concat
    df = pd.concat([df, minmax_df], axis=1)
    log_dataframe(df, log_message="  - Input SGEs:", indentation="  ")

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
    df_subset["x_offset_in_um"] = df_subset["x_offset"] / args.units_per_um if args.convert_to_um else np.nan
    df_subset["y_offset_in_um"] = df_subset["y_offset"] / args.units_per_um if args.convert_to_um else np.nan

    df_subset["new_xmin"] = df_subset["xmin"] + df_subset["x_offset"]
    df_subset["new_xmax"] = df_subset["xmax"] + df_subset["x_offset"]
    df_subset["new_ymin"] = df_subset["ymin"] + df_subset["y_offset"]
    df_subset["new_ymax"] = df_subset["ymax"] + df_subset["y_offset"]

    # * Convert to um (to allow offset_in_um in out_offsets file)
    df_subset["new_xmin_in_um"] = df_subset["new_xmin"] / args.units_per_um if args.convert_to_um else np.nan
    df_subset["new_xmax_in_um"] = df_subset["new_xmax"] / args.units_per_um if args.convert_to_um else np.nan
    df_subset["new_ymin_in_um"] = df_subset["new_ymin"] / args.units_per_um if args.convert_to_um else np.nan
    df_subset["new_ymax_in_um"] = df_subset["new_ymax"] / args.units_per_um if args.convert_to_um else np.nan

    log_dataframe(df_subset,  log_message="  - Offsets and new minmax per tile:", indentation="  ")

    # * save subset offsets & minmax
    out_subset = os.path.join(args.out_dir, args.out_subset)
    df_subset.to_csv(out_subset, sep="\t", index=False)
    logger.info(f"  - Offsets and new minmax per tile written to {out_subset}")

    # 4. Minmax for the combined SGEs
    logger.info(f"  "+"-"*20)
    logger.info(f"  - Calculating minmax for the combined SGEs...")

    # * in unit
    df_minmax_combined = {
        "xmin": 0,
        "xmax": col2x.sum(), 
        "ymin": 0,
        "ymax": row2y.sum()
    }
    df_minmax_combined = {k: round(v, args.precision) for k, v in df_minmax_combined.items()} 
    logger.info(f"  - Minmax for the combined SGEs (in unit): {df_minmax_combined}")

    # * in um
    if args.convert_to_um:
        df_minmax_combined["xmin"] = df_minmax_combined["xmin"] / args.units_per_um
        df_minmax_combined["xmax"] = df_minmax_combined["xmax"] / args.units_per_um
        df_minmax_combined["ymin"] = df_minmax_combined["ymin"] / args.units_per_um
        df_minmax_combined["ymax"] = df_minmax_combined["ymax"] / args.units_per_um
        df_minmax_combined = {k: round(v, args.precision) for k, v in df_minmax_combined.items()} 
        logger.info(f"  - Minmax for the combined SGEs (in um): {df_minmax_combined}")

    # * Write minmax for combined SGE
    out_minmax = os.path.join(args.out_dir, args.out_minmax)
    with open(out_minmax, "w") as f:
        f.writelines([f"{key}\t{value}\n" for key, value in df_minmax_combined.items()])
    
    logger.info(f"  - Minmax file written to {out_minmax}")

    # 5. feature
    logger.info(f"  "+"-"*20)
    logger.info(f"  - Combining feature counts...")

    def combine_ftr_across_sge(file_list, count_columns, feature_name_col, feature_id_col=None, debug=False):
        """Combines feature counts across multiple files."""
        df_ftr = pd.concat([pd.read_csv(file, sep="\t") for file in file_list])
        group_cols = [feature_name_col] + ([feature_id_col] if feature_id_col else [])
        df_ftr = df_ftr.groupby(group_cols)[count_columns].sum().reset_index()
        if debug:
            df_ftr[df_ftr[count_columns].sum(axis=1) == 0].to_csv("empty_feature.tsv", sep="\t", index=False)
        # Drop rows with all zero counts
        df_ftr = df_ftr[df_ftr[count_columns].sum(axis=1) > 0]
        return df_ftr

    out_ftr = os.path.join(args.out_dir, args.out_feature)
    in_ftrs = df["feature_path"].tolist()
    df_ftr_combined = combine_ftr_across_sge(in_ftrs, args.colnames_count, args.colname_feature_name, args.colname_feature_id, debug=args.debug)
    df_ftr_combined.to_csv(out_ftr, sep="\t", index=False, compression="gzip")
    logger.info(f"  - Feature file written to {out_ftr}")

    # 6. transcript
    logger.info(f"  "+"-"*20)
    logger.info(f"  - Combining transcripts...")

    def combine_transcript_across_sge(df, output_file, col_x, col_y, count_cols, out_cols, units_per_um, convert_to_um, precision, chunk_size=100000, debug=False):
        """Combines transcript data across multiple input files."""
        with gzip.open(output_file, "wt") as out_file:
            for idx, row in df.iterrows():
                try:
                    x_offset, y_offset = row["x_offset"], row["y_offset"]
                    transcript_path = row["transcript_path"]

                    logger.info(f" - Processing: {transcript_path} (Row: {row['row']}, Col: {row['col']})")
                    with gzip.open(transcript_path, "rt") as f:
                        for chunk in pd.read_csv(f, sep="\t", chunksize=chunk_size):
                            chunk[col_x] += x_offset
                            chunk[col_y] += y_offset
                
                            if convert_to_um:
                                chunk[[col_x, col_y]] = (chunk[[col_x, col_y]] / units_per_um).round(precision)
                            
                            chunk["index"] = idx
                            chunk = chunk[chunk[count_cols].sum(axis=1) > 0]  # Drop zero-count rows
                            assert all(col in chunk.columns for col in out_cols), f"Output columns not found in the chunk(cols: {chunk.columns}) "
                            chunk[out_cols].to_csv(out_file, sep="\t", index=False, header=(idx == 0))
                    if debug:
                        logger.info("Test mode enabled, stopping after first chunk.")
                        return
                except Exception as e:
                    raise ValueError(f"Error processing {transcript_path}: {e}")

    out_transcript = os.path.join(args.out_dir, args.out_transcript)
    combine_transcript_across_sge(df, out_transcript, 
                                args.colname_x, args.colname_y, args.colnames_count, out_cols=out_cols, 
                                units_per_um=args.units_per_um, convert_to_um=args.convert_to_um, precision=args.precision,
                                debug=args.debug)
    logger.info(f"Transcript file written to {out_transcript}")

    if args.debug:
        log_dataframe(df, log_message="Details:", indentation="  ")


if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
