import sys, os, gzip, argparse, logging, warnings, shutil, re, copy, time, pickle, inspect, warnings, json, yaml
import pandas as pd


def parse_arguments(_args):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}", 
                                    description=""" Generate a feature file for the SGE. """)
    parser.add_argument("--in-transcript", type=str, help="Input file.", required=True)
    parser.add_argument("--add-feature", action='store_true', help="Add feature to the input file.")
    parser.add_argument("--add-minmax", action='store_true', help="Add minmax to the input file.")

    # feature args
    parser.add_argument("--out-feature", type=str, help="Output file for feature.", default=None)
    parser.add_argument("--index-col", type=str, nargs='*', help="Column names to be used as index.", default=['gene_id', 'gene'])
    parser.add_argument("--count-col", type=str, nargs='*', help="Column names to be used as count.", default=['gn', 'gt', 'spl', 'unspl', 'ambig'])
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
    # # Determine the output file path if not provided
    # if args.out_feature is None:
    #     in_dir = os.path.dirname(args.in_transcript)
    #     in_id = os.path.basename(args.in_transcript).replace('.tsv.gz', '').replace(".transcripts", "").replace(".transcript", "")
    #     args.out_feature = os.path.join(in_dir, f"{in_id}.feature.tsv.gz")

    # Read the input transcript file
    with gzip.open(args.in_transcript, 'rt') as f:
        transcripts = pd.read_csv(f, sep='\t')

    # Ensure required columns are in the input
    missing_cols = set(args.index_col + args.count_col) - set(transcripts.columns)
    if missing_cols:
        raise ValueError(f"The following required columns are missing in the input: {', '.join(missing_cols)}")

    # Aggregate all feature columns in a single groupby operation
    feature_summary = transcripts.groupby(args.index_col, as_index=False)[args.count_col].sum()

    # Save the summarized DataFrame to a compressed TSV file
    with gzip.open(args.out_feature, 'wt') as f:
        feature_summary.to_csv(f, sep='\t', index=False)

    print(f"Feature saved to {args.out_feature}")

def sge_add_minmax(args):
    """Extract, scale, and compute min/max for X and Y coordinates."""
    # Determine the output file path if not provided
    assert args.out_minmax is not None, "When --add-minmax, --out-minmax must be provided."
    # if args.out_minmax is None:
    #     #raise ValueError("Output file for minmax must be provided.")
    #     in_dir = os.path.dirname(args.in_transcript)
    #     in_id = os.path.basename(args.in_transcript).replace('.tsv.gz', '').replace(".transcripts", "").replace(".transcript", "")
    #     args.out_minmax = os.path.join(in_dir, f"{in_id}.minmax.tsv")

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