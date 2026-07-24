import sys, os, gzip, argparse, inspect

from cartloader.utils.feature_filter import FeatureFilter


def flexopen(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Filter a molecule/transcript TSV by feature name, writing the surviving rows to a "
                    "new TSV. A standalone, on-demand version of the feature filtering that run_together "
                    "applies during a pipeline: use it to prepare a custom-filtered transcript file "
                    "without rerunning ingest.")

    inout = parser.add_argument_group("Input/Output")
    inout.add_argument('--in', dest='in_tsv', required=True, type=str,
                       help='Input transcript TSV (.tsv or .tsv.gz), with a header row')
    inout.add_argument('--out', required=True, type=str,
                       help='Output transcript TSV (.tsv or .tsv.gz); gzip-compressed if the name ends in .gz')

    col = parser.add_argument_group("Feature Column")
    col.add_argument('--colname-feature', type=str, default='gene',
                     help='Name of the feature/gene column in the input header (default: gene)')
    col.add_argument('--colidx-feature', type=int, default=None,
                     help='1-based index of the feature column. Overrides --colname-feature; use it when '
                          'the header does not name the column as expected.')

    ftr = parser.add_argument_group("Feature Filtering")
    FeatureFilter.add_arguments(ftr, "the output")

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def filter_molecules(_args):
    args = parse_arguments(_args)

    if not os.path.exists(args.in_tsv):
        sys.exit(f"ERROR: input file not found: {args.in_tsv} (--in)")
    if os.path.abspath(args.in_tsv) == os.path.abspath(args.out):
        sys.exit("ERROR: --in and --out must differ (filtering does not happen in place).")
    for flag, path in (("--include-feature-list", args.include_feature_list),
                       ("--exclude-feature-list", args.exclude_feature_list)):
        if path is not None and not os.path.exists(path):
            sys.exit(f"ERROR: file not found: {path} ({flag})")

    ftr_filter = FeatureFilter.from_args(args)
    if not ftr_filter.active:
        sys.exit("ERROR: no feature filter given. Provide at least one of --exclude-feature-regex, "
                 "--exclude-feature-list, --include-feature-regex, --include-feature-list.")

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    n_in = n_out = 0
    with flexopen(args.in_tsv, "rt") as rf, flexopen(args.out, "wt") as wf:
        header = rf.readline()
        if not header:
            sys.exit(f"ERROR: {args.in_tsv} is empty.")
        cols = header.rstrip("\n").split("\t")

        if args.colidx_feature is not None:
            if not (1 <= args.colidx_feature <= len(cols)):
                sys.exit(f"ERROR: --colidx-feature {args.colidx_feature} is out of range for the {len(cols)} "
                         f"columns in {args.in_tsv}.")
            icol = args.colidx_feature - 1
        else:
            if args.colname_feature not in cols:
                sys.exit(f"ERROR: feature column '{args.colname_feature}' not found in the header of "
                         f"{args.in_tsv}: {cols}. Set --colname-feature or --colidx-feature.")
            icol = cols.index(args.colname_feature)

        wf.write(header)
        for line in rf:
            toks = line.rstrip("\n").split("\t")
            if len(toks) <= icol:
                continue
            n_in += 1
            if ftr_filter.keep(toks[icol]):
                wf.write(line)
                n_out += 1

    print(f"Filtered {args.in_tsv}: {n_in} molecules read, {ftr_filter.n_dropped} dropped, "
          f"{n_out} written -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    filter_molecules(sys.argv[1:])
