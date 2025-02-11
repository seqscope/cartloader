import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd

from cartloader.utils.utils import log_dataframe

def parse_minmax(file_path):
    with open(file_path, "r") as f:
        return {k: float(v) for k, v in (line.strip().split("\t") for line in f)}

def combine_ftr_across_sge(file_list, colnames_count, colname_feature_name, colname_feature_id=None, test_mode=False):
    data_frames = [pd.read_csv(file, sep="\t") for file in file_list]
    df_ftr_concat = pd.concat(data_frames)
    group_columns = [colname_feature_name]
    if colname_feature_id and colname_feature_id in df_ftr_concat.columns:
        group_columns.append(colname_feature_id)
    df_ftr_combined = df_ftr_concat.groupby(group_columns)[colnames_count].sum().reset_index()
    # drop rows with all zero counts
    if test_mode:
        df_ftr_combined_empty= df_ftr_combined.loc[(df_ftr_combined[colnames_count] == 0).all(axis=1)]
        df_ftr_combined_empty.to_csv("empty_feature.tsv", sep="\t", index=False)
    df_ftr_combined = df_ftr_combined.loc[~(df_ftr_combined[colnames_count] == 0).all(axis=1)]
    return df_ftr_combined

def combine_transcript_across_sge(df, out_transcript, colx, coly, colnames_count, outcols, units_per_um, out_in_um, chunksize=100000, test_mode=False):
    logging.info("Merging transcripts...")
    with gzip.open(out_transcript, "wt") as out_file:
        for idx, row in df.iterrows():
            try:
                x_offset, y_offset = row["x_offset"], row["y_offset"]

                logging.info(f" - File: {row['transcript_path']}")
                logging.info(f"   row,col: {row['row']},{row['col']}")

                if out_in_um:
                    logging.info(f"   x y offsets (in um): {x_offset / units_per_um:.2f},{y_offset / units_per_um:.2f}")
                    logging.info(
                        f"   subset minmax(xmin,xmax,ymin,ymax) (in um): "
                        f"{row['xmin'] / units_per_um:.2f},{row['xmax'] / units_per_um:.2f},"
                        f"{row['ymin'] / units_per_um:.2f},{row['ymax'] / units_per_um:.2f}"
                    )
                else:
                    logging.info(f"   x y offsets (in unit): {x_offset},{y_offset}")
                    logging.info(
                        f"   subset minmax(xmin,xmax,ymin,ymax) (in unit): "
                        f"{row['xmin']},{row['xmax']},{row['ymin']},{row['ymax']}"
                    )

                with gzip.open(row["transcript_path"], "rt") as f:
                    for chunk in pd.read_csv(f, sep="\t", chunksize=chunksize):
                        chunk[colx] += x_offset
                        chunk[coly] += y_offset

                        if out_in_um:
                            chunk[[colx, coly]] = (chunk[[colx, coly]] / units_per_um).round(2)

                        chunk["index"] = idx
                        header = idx == 0 and chunk.index[0] == 0
                        # drop rows with all zero counts
                        chunk = chunk.loc[~(chunk[colnames_count] == 0).all(axis=1)]
                        chunk[outcols].to_csv(out_file, sep="\t", index=False, header=header)

                        if test_mode:
                            logging.info("Test mode is True. Exiting after processing the first chunk.")
                            return
            except Exception as e:
                logging.error(f"Error processing {row['transcript_path']}: {e}")
                raise

# def combine_transcript_across_sge(df, out_transcript, colx, coly, outcols, units_per_um, out_in_um, chunksize=100000, test_mode=False):
#     logging.info("Merging transcripts...")
#     with gzip.open(out_transcript, "wt") as out_file:
#         for idx, row in df.iterrows():
#             try:
#                 x_offset = row["x_offset"]
#                 y_offset = row["y_offset"]
#                 logging.info(f" - File: {row['transcript_path']}")
#                 logging.info(f"   row,col: {row['row']},{row['col']}") 
#                 if out_in_um:
#                     logging.info(f"   x y offsets (in um): {round(x_offset/units_per_um, 2)},{round(y_offset/units_per_um,2)}")
#                     logging.info(f"   subset minmax(xmin,xmax,ymin,ymax) (in um): {round(row['xmin']/units_per_um, 2)},{round(row['xmax']/units_per_um,2)},{round(row['ymin']/units_per_um,2)},{round(row['ymax']/units_per_um,2)}")
#                 else:
#                     logging.info(f"   x y offsets (in unit): {x_offset},{y_offset}")
#                     logging.info(f"   subset minmax(xmin,xmax,ymin,ymax) (in unit): {row['xmin']},{row['xmax']},{row['ymin']},{row['ymax']}")
#                 with gzip.open(row["transcript_path"], "rt") as f:
#                     for chunk in pd.read_csv(f, sep="\t", chunksize=chunksize):
#                         chunk[colx] += x_offset
#                         chunk[coly] += y_offset
#                         chunk["index"] = idx
#                         # convert into um with two decimal points if needed
#                         if out_in_um:
#                             chunk[colx] = chunk[colx] / units_per_um
#                             chunk[coly] = chunk[coly] / units_per_um
#                             chunk[colx] = chunk[colx].round(2)
#                             chunk[coly] = chunk[coly].round(2)
#                         header = idx == 0 and chunk.index[0] == 0
#                         chunk = chunk[outcols]
#                         chunk.to_csv(out_file, sep="\t", index=False, header=header)                        
#                         # If test mode is enabled, process only the first chunk
#                         if test_mode:
#                             logging.info("Test mode is True. Exiting after processing the first chunk.")
#                             break
#             except Exception as e:
#                 logging.error(f"Error processing {row['transcript_path']}: {e}")
#                 raise


def combine_sges_by_layout(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Combine SGE data by layout"
    )
    parser.add_argument("--in-tiles", type=str, nargs='*', default=[], help="List of input information in a specific format.")
    parser.add_argument("--out-dir", type=str, help="Output directory.")
    #parser.add_argument("--outid", type=str, default=None, help="Output ID.")
    parser.add_argument('--out-transcript', type=str, default="transcripts.unsorted.tsv.gz", help='Output for SGE combination. The compressed transcript-indexed SGE file in TSV format (default: transcripts.unsorted.tsv.gz).')
    parser.add_argument('--out-minmax', type=str, default="coordinate_minmax.tsv", help='Output for SGE combination. The coordinate minmax TSV file (default: coordinate_minmax.tsv).')
    parser.add_argument('--out-feature', type=str, default="feature.clean.tsv.gz", help='Output for SGE combination. The compressed UMI count per gene TSV file (default: feature.clean.tsv.gz).')
    parser.add_argument("--colnames-count", type=str, nargs='*', help="Columns to sum (default: count).", default=['count'])
    parser.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    parser.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    parser.add_argument('--colname-x', type=str, default="X", help='X column name (default: X)')
    parser.add_argument('--colname-y', type=str, default="Y", help='Y column name (default: Y)')
    parser.add_argument('--units-per-um', type=float, default=1.0, help='Units per um in the input transcript tsv files (default: 1.0)')
    parser.add_argument('--minmax-in-um', action='store_true', help='Input minmax is in um while input transcript is based on unit.')
    parser.add_argument('--convert-to-um', action='store_true', help='Convert output to um.')
    parser.add_argument('--test-mode', action='store_true', help='Test mode.')
    # parser.add_argument("--unit-ids", type=str, nargs='*', default=[], help="(For test) List of input unit ids.")
    # parser.add_argument("--analy-dir", type=str, help="(For test) Analysis directory.")
    args = parser.parse_args(_args)

    os.makedirs(args.out_dir, exist_ok=True) 
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        #filename=os.path.join(args.out_dir, f"{args.outid}.log" if args.outid else "sge_combine.log"),
        filename=os.path.join(args.out_dir, "sge_stitch.log"),
        filemode='w'
    )

    input_list = args.in_tiles

    # log
    # logging.info("="*30)
    # logging.info(f"Commands:")
    # cmd=" ".join([f"cartloader combine_sges_by_layout",
    #     f"--input {' '.join(input_list)} ",
    #     f"--out-dir {args.out_dir}",
    #     f"--outid {args.outid} " if args.outid else "",
    #     f"--colnames-count {' '.join(args.colnames_count)} ",
    #     f"--colname-feature-name {args.colname_feature_name} ",
    #     f"--colname-feature-id {args.colname_feature_id} " if args.colname_feature_id else "",
    #     f"--colname-x {args.colname_x} ",
    #     f"--colname-y {args.colname_y} ",
    #     f"--units-per-um {args.units_per_um} ",
    #     '--minmax-in-um' if args.minmax_in_um else '',
    #     '--convert-to-um' if args.convert_to_um else '',
    #     f'--test-mode' if args.test_mode else ''])
    # logging.info(cmd)
    # logging.info("="*30)

    logging.info(f"  - Output directory: {args.out_dir}")

    outcols=[args.colname_x, args.colname_y]
    outcols.extend([args.colname_feature_name]) if args.colname_feature_id is None else outcols.extend([args.colname_feature_name, args.colname_feature_id])
    outcols.extend(args.colnames_count)

    df = pd.DataFrame(input_list, columns=["input"])
    df = df["input"].str.split(",", expand=True)
    df.columns = ["transcript_path", "feature_path", "minmax_path", "row", "col"]

    log_dataframe(df, log_message="Input SGEs:", indentation="  ")

    df["row"] = df["row"].astype(int)
    df["col"] = df["col"].astype(int)

    minmax_data = df["minmax_path"].apply(parse_minmax)
    minmax_df = pd.DataFrame(minmax_data.tolist())

    if args.minmax_in_um and args.units_per_um != 1.0:
        minmax_df["xmin"] = minmax_df["xmin"] * args.units_per_um
        minmax_df["xmax"] = minmax_df["xmax"] * args.units_per_um
        minmax_df["ymin"] = minmax_df["ymin"] * args.units_per_um
        minmax_df["ymax"] = minmax_df["ymax"] * args.units_per_um

    df = pd.concat([df, minmax_df], axis=1)

    df["x"] = df["xmax"] - df["xmin"]
    df["y"] = df["ymax"] - df["ymin"]

    row2y = df.groupby("row")["y"].max()
    col2x = df.groupby("col")["x"].max()

    col_offsets = col2x.cumsum() - col2x / 2
    row_offsets = row2y.cumsum() - row2y / 2

    df["x_offset"] = df["col"].map(col_offsets) - (df["xmin"] + df["xmax"]) / 2
    df["y_offset"] = df["row"].map(row_offsets) - (df["ymin"] + df["ymax"]) / 2

    df["x_min_subset"] = df["xmin"] + df["x_offset"]
    df["x_max_subset"] = df["xmax"] + df["x_offset"]
    df["y_min_subset"] = df["ymin"] + df["y_offset"]
    df["y_max_subset"] = df["ymax"] + df["y_offset"]

    # Minmax
    # - in unit
    df_minmax_combined = {
        "xmin": 0,
        "xmax": col2x.sum(),
        "ymin": 0,
        "ymax": row2y.sum()
    }
    # - in um
    if args.convert_to_um:
        df_minmax_combined["xmin"] = df_minmax_combined["xmin"] / args.units_per_um
        df_minmax_combined["xmax"] = df_minmax_combined["xmax"] / args.units_per_um
        df_minmax_combined["ymin"] = df_minmax_combined["ymin"] / args.units_per_um
        df_minmax_combined["ymax"] = df_minmax_combined["ymax"] / args.units_per_um
    
    df_minmax_combined = {k: round(v, 2) for k, v in df_minmax_combined.items()}

    #out_minmax = os.path.join(args.out_dir, f"{args.outid}.coordinate_minmax.tsv" if args.outid else "coordinate_minmax.tsv")
    out_minmax = os.path.join(args.out_dir, args.out_minmax)
    with open(out_minmax, "w") as f:
        f.writelines([f"{key}\t{value}\n" for key, value in df_minmax_combined.items()])
    logging.info(f"Minmax file written to {out_minmax}")

    # feature
    #out_ftr = os.path.join(args.out_dir, f"{args.outid}.feature.tsv.gz" if args.outid else "feature.clean.tsv.gz")
    out_ftr = os.path.join(args.out_dir, args.out_feature)
    in_ftrs = df["feature_path"].tolist()
    df_ftr_combined = combine_ftr_across_sge(in_ftrs, args.colnames_count, args.colname_feature_name, args.colname_feature_id, test_mode=args.test_mode)
    df_ftr_combined.to_csv(out_ftr, sep="\t", index=False, compression="gzip")
    logging.info(f"Feature file written to {out_ftr}")

    #out_transcript = os.path.join(args.out_dir, f"{args.outid}.transcripts.unsorted.tsv.gz" if args.outid else "transcripts.unsorted.tsv.gz")
    out_transcript = os.path.join(args.out_dir, args.out_transcript)
    combine_transcript_across_sge(df, out_transcript, args.colname_x, args.colname_y, args.colnames_count, outcols=outcols, units_per_um=args.units_per_um, out_in_um=args.convert_to_um, test_mode=args.test_mode)
    logging.info(f"Transcript file written to {out_transcript}")
    
    if args.test_mode:
        log_dataframe(df, log_message="Details:", indentation="  ")


if __name__ == "__main__":
    func_name = os.path.splitext(os.path.basename(__file__))[0]
    func = getattr(sys.modules[__name__], func_name)
    func(sys.argv[1:])
