import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
import subprocess
import numpy as np

from cartloader.utils.utils import log_dataframe, create_custom_logger

def sge_combine_tiles(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Combine SGE data by layout. Note the transcript/feature/minmax files should be in the same resolution.",
    )
    parser.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of input information in a specific format: <transcript_path>,<feature_path>,<minmax_path>,<row>,<col>.")
    parser.add_argument('--in-tile-minmax', type=str, default=None, help="Path to the input offsets file, in which the following columns are required: row, col, x_offset_unit, y_offset_unit, global_xmin_um, global_xmax_um, global_ymin_um, global_ymax_um")
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
    parser.add_argument("--units-per-um", type=float, default=1.0, help="Define the scaling factor for unit conversion from coordinate to um (default: 1.0)")
    parser.add_argument("--precision", type=int, default=2, help="Precision for the output minmax and transcript files (default: 2)")
    parser.add_argument('--debug', action='store_true', help='Test mode.')
    parser.add_argument('--log', action='store_true', default=False, help='Write log to file')

    args = parser.parse_args(_args)

    # outdir
    if os.path.isfile(args.out_dir):
        raise NotADirectoryError(f"Output path exists as a file, not a directory: {args.out_dir}")

    os.makedirs(args.out_dir, exist_ok=True) 

    # log
    logger = create_custom_logger(__name__, os.path.join(args.out_dir, f"{args.out_transcript}.log") if args.log else None)
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

    # read in df
    df = pd.DataFrame(args.in_tiles, columns=["input"])
    df = df["input"].str.split(",", expand=True)
    df.columns = ["transcript_path", "feature_path", "minmax_path", "row", "col"]
    df["row"] = df["row"].astype(int)
    df["col"] = df["col"].astype(int)
    log_dataframe(df, msg="  - Input SGEs:")
    df["layout_idx"] = df.apply(lambda x: f"{int(x['row'])},{int(x['col'])}", axis=1)

    # read in_offsets
    logger.info(f"  - Reading in offsets...")
    df_tile_coords = pd.read_csv(args.in_tile_minmax, sep="\t", header=0)
    assert df.shape[0] == df_tile_coords.shape[0], f"Found mismatch number of tiles between --in-tiles ({df.shape[0]}) and --in-tile-minmax ({df_tile_coords.shape[0]})."
    df_tile_coords["layout_idx"] = df_tile_coords.apply(lambda x: f"{int(x['row'])},{int(x['col'])}", axis=1)

    # update df: input with offsets
    df=df.merge(df_tile_coords[["layout_idx", "x_offset_unit", "y_offset_unit"]], on="layout_idx", how="left")
    
    # 3. Minmax for the combined SGEs
    logger.info(f"  - Get minmax for the combined SGEs ...")
    # make sure the row and col are 
    global_xmin_um = df_tile_coords["global_xmin_um"].min()
    global_xmax_um = df_tile_coords["global_xmax_um"].max()
    global_ymin_um = df_tile_coords["global_ymin_um"].min()
    global_ymax_um = df_tile_coords["global_ymax_um"].max()

    df_minmax_combined = {
        "xmin": global_xmin_um,
        "xmax": global_xmax_um,
        "ymin": global_ymin_um,
        "ymax": global_ymax_um
    }
    df_minmax_combined = {k: round(v, args.precision) for k, v in df_minmax_combined.items()} 
    logger.info(f"  - Minmax for the combined SGEs (in um): {df_minmax_combined}")

    out_minmax = os.path.join(args.out_dir, args.out_minmax)
    with open(out_minmax, "w") as f:
        f.writelines([f"{key}\t{value}\n" for key, value in df_minmax_combined.items()])
    
    logger.info(f"  - Minmax file written to {out_minmax}")

    # 4. feature
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

    # 5. transcript
    logger.info(f"  "+"-"*20)
    logger.info(f"  - Combining transcripts...")
    logger.info(f"  - Units per um: {args.units_per_um}")
    logger.info(f"  - Precision: {args.precision}")

    def combine_transcript_across_sge(df, output_file, col_x, col_y, count_cols, out_cols, units_per_um, precision, chunk_size=100000, debug=False):
        """Combines transcript data across multiple input files."""
        if os.path.exists(output_file):
            logger.info(f"  - Output file already exists, removing: {output_file}")
            os.remove(output_file)
        with gzip.open(output_file, "wt") as out_file:
            chunk_gidx = 0 # global chunk index for handling the header
            for idx, row in df.iterrows():
                try:
                    x_offset, y_offset = row["x_offset_unit"], row["y_offset_unit"]
                    transcript_path = row["transcript_path"]
                    logger.info(f"  - Processing: {transcript_path} (Row: {row['row']}, Col: {row['col']}): x_offset={x_offset}, y_offset={y_offset}")
                    with gzip.open(transcript_path, "rt") as in_f:
                        for chunk in pd.read_csv(in_f, sep="\t", chunksize=chunk_size):
                            chunk[col_x] += float(x_offset)
                            chunk[col_y] += float(y_offset)
                            chunk[[col_x, col_y]] = (chunk[[col_x, col_y]] / units_per_um).round(precision)
                            chunk = chunk[chunk[count_cols].sum(axis=1) > 0]  # Drop zero-count rows
                            assert all(col in chunk.columns for col in out_cols), f"Output columns not found in the chunk(cols: {chunk.columns}) "
                            chunk[out_cols].to_csv(out_file, sep="\t", index=False, header=(chunk_gidx == 0))
                            chunk_gidx += 1
                    if debug:
                        logger.info("Test mode enabled, stopping after first chunk.")
                        return
                except Exception as e:
                    raise ValueError(f"Error processing {transcript_path}: {e}")

    out_transcript = os.path.join(args.out_dir, args.out_transcript)
    combine_transcript_across_sge(df, out_transcript, 
                                args.colname_x, args.colname_y, args.colnames_count, out_cols=out_cols, 
                                units_per_um=args.units_per_um,  precision=args.precision,
                                debug=args.debug)
    logger.info(f"Transcript file written to {out_transcript}")

    if args.debug:
        log_dataframe(df, msg="Details:", indentation="  ")


if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
