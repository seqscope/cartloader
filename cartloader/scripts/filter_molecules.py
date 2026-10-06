import sys, os, argparse, inspect, shlex, subprocess

from cartloader.utils.feature_filter import FeatureFilter


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

    env = parser.add_argument_group("ENV Parameters")
    env.add_argument('--gzip', type=str, default="gzip",
                     help='Program used to (de)compress .gz input/output, run as a separate process so it '
                          'does not block the Python loop. For a large file use "pigz" (or "pigz -p 8") for '
                          'multi-threaded compression (default: gzip).')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def _open_reader(path, gzip_prog):
    """(proc, binary stream) for reading `path`. A .gz input is decompressed by an external
    process so decompression overlaps with — and runs on a different core than — the filter
    loop; a plain file is opened directly."""
    if path.endswith(".gz"):
        try:
            proc = subprocess.Popen(shlex.split(gzip_prog) + ["-dc", path], stdout=subprocess.PIPE)
        except FileNotFoundError:
            sys.exit(f"ERROR: cannot run decompressor '{gzip_prog}' (--gzip) for {path}")
        return proc, proc.stdout
    return None, open(path, "rb")


def _open_writer(path, gzip_prog):
    """(proc, binary write stream, file) for writing `path`. A .gz output is piped through an
    external compressor (use pigz for multi-threaded compression); a plain file is written
    directly. `file` is the underlying .gz file handle to close after the compressor exits."""
    if path.endswith(".gz"):
        fout = open(path, "wb")
        try:
            proc = subprocess.Popen(shlex.split(gzip_prog), stdin=subprocess.PIPE, stdout=fout)
        except FileNotFoundError:
            fout.close()
            sys.exit(f"ERROR: cannot run compressor '{gzip_prog}' (--gzip) for {path}")
        return proc, proc.stdin, fout
    return None, open(path, "wb"), None


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

    rproc, rf = _open_reader(args.in_tsv, args.gzip)
    wproc, wf, wfout = _open_writer(args.out, args.gzip)

    n_in = n_out = 0
    try:
        header = rf.readline()
        if not header:
            sys.exit(f"ERROR: {args.in_tsv} is empty.")
        cols = header.rstrip(b"\r\n").split(b"\t")

        if args.colidx_feature is not None:
            if not (1 <= args.colidx_feature <= len(cols)):
                sys.exit(f"ERROR: --colidx-feature {args.colidx_feature} is out of range for the {len(cols)} "
                         f"columns in {args.in_tsv}.")
            icol = args.colidx_feature - 1
        else:
            try:
                icol = cols.index(args.colname_feature.encode())
            except ValueError:
                names = [c.decode("utf-8", "replace") for c in cols]
                sys.exit(f"ERROR: feature column '{args.colname_feature}' not found in the header of "
                         f"{args.in_tsv}: {names}. Set --colname-feature or --colidx-feature.")
        wf.write(header)

        # The kept/dropped decision depends only on the feature name, whose cardinality is
        # small (the gene panel), so memoize it: the regex/set logic in FeatureFilter.keep
        # runs once per distinct gene, and the per-molecule hot path is a dict lookup plus a
        # raw-bytes write. Work stays in bytes (no decode of the whole row); only the feature
        # token, and only on a cache miss, is decoded. (This makes FeatureFilter.n_dropped
        # count distinct dropped genes, not rows, so molecule totals are derived from n_out.)
        keep = ftr_filter.keep
        write = wf.write
        cache = {}
        nsplit = icol + 1
        for line in rf:
            parts = line.split(b"\t", nsplit)
            if len(parts) <= icol:
                continue
            n_in += 1
            feat = parts[icol]
            if len(parts) == icol + 1:        # feature is the last field -> drop its newline
                feat = feat.rstrip(b"\r\n")
            decision = cache.get(feat)
            if decision is None:
                decision = keep(feat.decode("utf-8", "replace"))
                cache[feat] = decision
            if decision:
                write(line)
                n_out += 1
    finally:
        # Close the write side first (EOF to the compressor), then the read side, and surface
        # a non-zero exit from either helper process.
        wf.close()
        w_rc = wproc.wait() if wproc is not None else 0
        if wfout is not None:
            wfout.close()
        rf.close()
        r_rc = rproc.wait() if rproc is not None else 0
    if w_rc != 0:
        sys.exit(f"ERROR: compressor '{args.gzip}' failed (exit {w_rc}) writing {args.out}")
    if r_rc != 0:
        sys.exit(f"ERROR: decompressor '{args.gzip}' failed (exit {r_rc}) reading {args.in_tsv}")

    print(f"Filtered {args.in_tsv}: {n_in} molecules read, {n_in - n_out} dropped, "
          f"{n_out} written -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    filter_molecules(sys.argv[1:])
