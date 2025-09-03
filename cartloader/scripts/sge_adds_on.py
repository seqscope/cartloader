import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd
from collections import defaultdict
import math

"""
cartloader sge_adds_on \
    --in-transcript /nfs/turbo/umms-leeju/nova/v2/analysis/n14-hm2tk-t07a-mouse-1bbd0/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207/preprocess/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207.transcripts.tsv.gz \
    --out-feature /nfs/turbo/umms-leeju/nova/v2/analysis/n14-hm2tk-t07a-mouse-1bbd0/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207/preprocess/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207.feature2.tsv.gz \
    --colname-feature-id gene_id --add-feature

cartloader sge_adds_on \
    --in-transcript /nfs/turbo/umms-leeju/nova/v2/analysis/n14-hm2tk-t07a-mouse-1bbd0/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207/preprocess/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207.transcripts.tsv.gz \
    --out-minmax /nfs/turbo/umms-leeju/nova/v2/analysis/n14-hm2tk-t07a-mouse-1bbd0/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207/preprocess/n14-hm2tk-t07a-mouse-1bbd0-mask-a3207.gn.raw.coordinate_minmax2.tsv\
    --add-minmax --mu-scale 1000
"""

def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description="""Generate a feature or a minmax file for the SGE in FICTURE-compatible format. """)
    parser.add_argument("--in-transcript", type=str, help="Input file.", required=True)
    parser.add_argument("--chunk-size", type=int, default=10**6, help="Chunk size for processing (default: 10^6).")
    # feature args
    parser.add_argument("--add-feature", action='store_true', help="Create a feature file based on the input file.")
    parser.add_argument("--out-feature", type=str, default=None, help="Path to output feature file.")
    parser.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    parser.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    parser.add_argument("--colname-count", type=str, default="gn", help="Comma-separated column names for count (default: gn)")
    # minmax args
    parser.add_argument("--add-minmax", action='store_true', help="Create a minmax file based on the input file.")
    parser.add_argument("--out-minmax", type=str, default=None, help="Path to output minmax file.")
    parser.add_argument("--mu-scale", type=float, default=1, help="Scale factor for X and Y coordinates (default: 1).")
    parser.add_argument("--colname-x", type=str, default="X", help="X coordinate column name (default: X).")
    parser.add_argument("--colname-y", type=str, default="Y", help="Y coordinate column name (default: Y).")
    parser.add_argument("--minmax-format", type=str, default="col", choices=["row", "col"], help="Format of the output minmax file. Options: 'row' outputs a single-line table with xmin, xmax, ymin, ymax as column headers; 'col' outputs a two-column format with minmax names and their values listed vertically (default: col).")
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)

def sge_add_feature_by_chunk(args):
    """Generate a summarized feature file using chunked processing to avoid OOM."""

    assert args.out_feature is not None, "When --add-feature, --out-feature must be provided."
    
    ftr_cols = [args.colname_feature_name, args.colname_feature_id] if args.colname_feature_id else [args.colname_feature_name]
    count_cols = [args.colname_count]

    aggregation = defaultdict(lambda: [0] * len(count_cols))  

    chunk_idx = 0
    with gzip.open(args.in_transcript, 'rt') as f:
        reader = pd.read_csv(f, sep='\t', chunksize=args.chunk_size)
        for chunk in reader:
            # Check for missing columns in first chunk
            if chunk_idx == 0:
                missing_cols = set(ftr_cols + count_cols) - set(chunk.columns)
                if missing_cols:
                    raise ValueError(f"The following required columns are missing in the input: {', '.join(missing_cols)}")
                chunk_idx = chunk_idx + 1 

            grouped = chunk.groupby(ftr_cols)[count_cols].sum()
            for idx, row in grouped.iterrows():
                if not isinstance(idx, tuple):
                    idx = (idx,)  # make single element group keys into tuple
                for i, val in enumerate(row):
                    aggregation[idx][i] += val

    # Create the final DataFrame
    agg_data = [list(key) + vals for key, vals in aggregation.items()]
    feature_summary = pd.DataFrame(agg_data, columns=ftr_cols + count_cols)

    # Save the summarized DataFrame to a compressed TSV file
    with gzip.open(args.out_feature, 'wt') as f:
        feature_summary.to_csv(f, sep='\t', index=False)

    print(f"Feature saved to {args.out_feature}")


def sge_add_minmax_by_chunk(args):
    """Extract, scale, and compute min/max for X and Y coordinates using chunked processing."""
    assert args.out_minmax is not None, "When --add-minmax, --out-minmax must be provided."

    xmin = ymin = math.inf
    xmax = ymax = -math.inf

    chunk_idx = 0
    with gzip.open(args.in_transcript, 'rt') as f:
        reader = pd.read_csv(f, sep="\t", chunksize=args.chunk_size)
        for chunk in reader:
            if chunk_idx == 0:
                # Check for X and Y columns
                assert args.colname_x in chunk.columns, f"Input file does not contain the required X column: {args.colname_x}."
                assert args.colname_y in chunk.columns, f"Input file does not contain the required Y column: {args.colname_y}."

            # drop the transcript with 0 count
            chunk = chunk[chunk[args.colname_count] != 0]

            x_scaled = chunk[args.colname_x] / args.mu_scale
            y_scaled = chunk[args.colname_y] / args.mu_scale

            xmin = min(xmin, x_scaled.min())
            xmax = max(xmax, x_scaled.max())
            ymin = min(ymin, y_scaled.min())
            ymax = max(ymax, y_scaled.max())

    # Prepare output
    minmax_dict = {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax}
    if args.minmax_format == "row":
        minmax_df = pd.DataFrame(minmax_dict, index=[0])
        with gzip.open(args.out_minmax, "wt") as f:
            minmax_df.to_csv(f, sep="\t", index=False)
    elif args.minmax_format == "col":
        with open(args.out_minmax, "w") as f:
            for key, value in minmax_dict.items():
                f.write(f"{key}\t{value}\n")

    print(f"Minmax saved to {args.out_minmax}")

def sge_adds_on(_args):
    """Generate a feature file for the SGE."""
    args = parse_arguments(_args)

    assert args.add_feature or args.add_minmax, "Either --add-feature or --add-minmax must be enabled."

    if args.add_feature:
        sge_add_feature_by_chunk(args)
    if args.add_minmax:
        sge_add_minmax_by_chunk(args)


if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])