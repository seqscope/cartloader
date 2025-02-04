import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd


def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description=""" Generate a feature or a minmax file for the SGE in FICTURE-compatible format. """)
    parser.add_argument("--in-transcript", type=str, help="Input file.", required=True)
    parser.add_argument("--add-feature", action='store_true', help="Add feature to the input file.")
    parser.add_argument("--add-minmax", action='store_true', help="Add minmax to the input file.")

    # feature args
    parser.add_argument("--out-feature", type=str, help="Output file for feature.", default=None)
    parser.add_argument('--colname-feature-name', type=str, default='gene', help='Feature name column (default: gene)')
    parser.add_argument('--colname-feature-id', type=str, default=None, help='Feature ID column (default: None)')
    parser.add_argument("--colnames-count", type=str, default="gn,gt,spl,unspl,ambig", help="Comma-separate column names for count (default: gn,gt,spl,unspl,ambig)")

    # minmax args
    parser.add_argument("--out-minmax", type=str, help="Output file for minmax.", default=None)
    parser.add_argument("--mu-scale", type=float, help="Scale factor for X and Y coordinates.", default=1)
    parser.add_argument("--minmax-format", type=str, help="Type of minmax to compute.", default="col")

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def sge_add_feature(args):
    """Generate a summarized feature file."""
    assert args.out_feature is not None, "When --add-feature, --out-feature must be provided."

    with gzip.open(args.in_transcript, 'rt') as f:
        transcripts = pd.read_csv(f, sep='\t')

    ftr_cols=[args.colname_feature_name, args.colname_feature_id] if args.colname_feature_id else [args.colname_feature_name]
    count_cols = args.colnames_count.split(",")

    missing_cols = set(ftr_cols + count_cols) - set(transcripts.columns)
    if missing_cols:
        raise ValueError(f"The following required columns are missing in the input: {', '.join(missing_cols)}")

    # Aggregate all feature columns in a single groupby operation
    feature_summary = transcripts.groupby(ftr_cols, as_index=False)[count_cols].sum()

    # Save the summarized DataFrame to a compressed TSV file
    with gzip.open(args.out_feature, 'wt') as f:
        feature_summary.to_csv(f, sep='\t', index=False)

    print(f"Feature saved to {args.out_feature}")

def sge_add_minmax(args):
    """Extract, scale, and compute min/max for X and Y coordinates."""
    assert args.out_minmax is not None, "When --add-minmax, --out-minmax must be provided."

    with gzip.open(args.in_transcript, 'rt') as f:
        transcripts = pd.read_csv(f, sep="\t")

    if "X" not in transcripts.columns or "Y" not in transcripts.columns:
        raise ValueError("Input file must contain 'X' and 'Y' columns.")

    transcripts["X_scaled"] = transcripts["X"] / args.mu_scale
    transcripts["Y_scaled"] = transcripts["Y"] / args.mu_scale

    xmin = transcripts["X_scaled"].min()
    xmax = transcripts["X_scaled"].max()
    ymin = transcripts["Y_scaled"].min()
    ymax = transcripts["Y_scaled"].max()

    # Return results as a dictionary
    minmax_dict={"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax}
    if args.minmax_format == "row":
        minmax_dict = pd.DataFrame(minmax_dict, index=[0])
        with gzip.open(args.out_minmax, "wt") as f:
            minmax_dict.to_csv(f, sep="\t", index=False)
    elif args.minmax_format == "col":
        with open(args.out_minmax, "w") as f:
            for key, value in minmax_dict.items():
                f.write(f"{key}\t{value}\n")
    print(f"Minmax saved to {args.out_minmax}")

def sge_adds_on(_args):
    """Generate a feature file for the SGE."""
    args = parse_arguments(_args)
    if args.add_feature:
        sge_add_feature(args)
    if args.add_minmax:
        sge_add_minmax(args)


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