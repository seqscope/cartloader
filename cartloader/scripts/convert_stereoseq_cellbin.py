import sys, os, gzip, argparse, inspect

# Stereo-seq cell-bin GEM -> the pixel TSV that `spatula pixel2sptsv` consumes.
#
# `saw convert gef2gem --cellbin-gef ... --cellbin-gem ...` writes a text GEM whose
# rows are (geneID, geneName, x, y, MIDCount, ExonCount, CellID): the subset of MIDs
# that cell segmentation placed inside a cell. Its coordinates cannot be mapped back
# onto the bin1 GEM row-by-row, so this file is NOT merged into the pixel transcript;
# it is converted here into a standalone 5-column TSV
#
#     X <TAB> Y <TAB> gene <TAB> count <TAB> cell_id        (no header)
#
# which run_together hands to run_ficture2_multi_cells as a `--tsv-list` entry. That
# column order matches the tool's --colidx-* defaults (1,2,3,4,5).
#
# Counts come from `ExonCount` (the exonic subset of the MIDs), matching the bin1
# ingest; zero-count rows are dropped rather than carried through as empty entries.
#
# The cell clusters derived from this file are projected onto the pixel-level model
# trained from the bin1 transcript, so the feature names here MUST match the ones
# sge_convert wrote for bin1 (both default to the GEM's `geneName` column, i.e. the
# gene symbol). Pass --check-features to verify that against the pixel run's feature
# file instead of discovering the mismatch as a run of near-empty cells.


def flexopen(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Convert a Stereo-seq cell-bin GEM (from 'saw convert gef2gem --cellbin-gem') "
                    "into the headerless X/Y/gene/count/cell_id TSV used for cell-level analysis.")

    inout = parser.add_argument_group("Input/Output Parameters")
    inout.add_argument('--in-gem', required=True, type=str, help='Path to the input cell-bin GEM (.gem or .gem.gz)')
    inout.add_argument('--out', required=True, type=str, help='Path to the output TSV (headerless; X, Y, gene, count, cell_id)')

    key = parser.add_argument_group("Key Parameters")
    key.add_argument('--units-per-um', type=float, default=2.0,
                     help='Coordinate units per um in the GEM. Stereo-seq bin1 coordinates are in '
                          '0.5um units, so 2 converts them to um (default: 2.0). Must match the '
                          '--units-per-um used by sge_convert for the bin1 GEM.')
    key.add_argument('--precision-um', type=int, default=2, help='Digits kept for the output coordinates (default: 2)')
    key.add_argument('--check-features', type=str, default=None,
                     help="Path to the pixel run's feature file (sge_convert's feature.clean.tsv.gz). "
                          'The cell clusters are projected onto a model trained on those features, so '
                          'this run fails if fewer than --min-feature-overlap of them appear here.')
    key.add_argument('--min-feature-overlap', type=float, default=0.5,
                     help='Fraction of --check-features features that must be present (default: 0.5)')
    key.add_argument('--check-features-colname', type=str, default='gene',
                     help="Feature-name column in --check-features (default: gene). This is "
                          "sge_convert's OUTPUT column name (--colname-feature-name), not the GEM's "
                          "input column; falls back to the first column if absent.")
    key.add_argument('--allow-offset', action='store_true', default=False,
                     help='Proceed even when the GEM header declares a non-zero OffsetX/OffsetY. By '
                          'default a non-zero offset is an error: the offset would shift the cells '
                          'relative to the pixel transcript and the registered histology.')

    incol = parser.add_argument_group("Input Column Parameters")
    incol.add_argument('--colname-feature', type=str, default='geneName',
                       help='Feature column in --in-gem (default: geneName, the gene symbol; must match '
                            'the bin1 ingest). Symbols are what CartoScope displays and what the '
                            'exclude-feature regexes match, so prefer them over geneID (Ensembl ids).')
    incol.add_argument('--colname-x', type=str, default='x', help='X column in --in-gem (default: x)')
    incol.add_argument('--colname-y', type=str, default='y', help='Y column in --in-gem (default: y)')
    incol.add_argument('--colname-count', type=str, default='ExonCount',
                       help='Count column in --in-gem (default: ExonCount, the exonic subset of the '
                            'MIDs; must match the bin1 ingest). Rows whose count is zero are dropped.')
    incol.add_argument('--colname-cell', type=str, default='CellID', help='Cell-id column in --in-gem (default: CellID)')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)


def read_check_features(path, colname_feature):
    """Read the feature names from sge_convert's feature file (a TSV whose header
    names the feature column with sge_convert's OUTPUT name, 'gene' by default)."""
    if not os.path.exists(path):
        sys.exit(f"ERROR: --check-features file not found: {path}\n"
                 f"       This should be the feature file sge_convert wrote for the bin1 GEM "
                 f"(--out-feature, default 'feature.clean.tsv.gz'). Point --check-features at it, "
                 f"or drop the flag to skip the feature-naming check.")
    names = set()
    with flexopen(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = header.index(colname_feature) if colname_feature in header else 0
        for line in f:
            toks = line.rstrip("\n").split("\t")
            if len(toks) > idx:
                names.add(toks[idx])
    return names


def convert_stereoseq_cellbin(_args):
    args = parse_arguments(_args)

    fmt = f"%.{args.precision_um}f" if args.precision_um >= 0 else "%.2f"
    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    seen_features = set()
    n_in = n_out = n_zero = 0

    with flexopen(args.in_gem, "rt") as rf:
        # The '#Key=Value' preamble carries the coordinate origin. A non-zero offset
        # means the body coordinates are relative, so they would not line up with the
        # bin1 transcript (whose own offset sge_convert does not apply either).
        header = None
        for line in rf:
            if line.startswith("#"):
                key, _, val = line[1:].strip().partition("=")
                if key in ("OffsetX", "OffsetY") and val.strip() not in ("", "0"):
                    msg = (f"{args.in_gem} declares {key}={val.strip()}. Non-zero offsets are not "
                           f"applied here, so the cells would be shifted relative to the pixel "
                           f"transcript and the histology. Re-export without the offset, or pass "
                           f"--allow-offset if you know the offset is already accounted for.")
                    if not args.allow_offset:
                        sys.exit(f"ERROR: {msg}")
                    print(f"WARNING: {msg}", file=sys.stderr)
                continue
            header = line.rstrip("\n").split("\t")
            break
        if header is None:
            sys.exit(f"ERROR: {args.in_gem} has no column header line (only '#' comments).")

        col = {}
        for label, name in (("ftr", args.colname_feature), ("x", args.colname_x), ("y", args.colname_y),
                            ("cnt", args.colname_count), ("cell", args.colname_cell)):
            if name not in header:
                sys.exit(f"ERROR: column '{name}' not found in the header of {args.in_gem}: {header}")
            col[label] = header.index(name)
        ncol = max(col.values()) + 1

        with open(args.out, "wt") as wf:
            for line in rf:
                toks = line.rstrip("\n").split("\t")
                if len(toks) < ncol:
                    continue
                n_in += 1
                # Drop zero-count rows. With ExonCount as the count column a large
                # share of rows are 0 (a MID with no exonic overlap); carrying them
                # through would add cells and features made entirely of empty counts.
                cnt = toks[col["cnt"]]
                if not float(cnt) > 0:
                    n_zero += 1
                    continue
                gene = toks[col["ftr"]]
                seen_features.add(gene)
                x = fmt % (float(toks[col["x"]]) / args.units_per_um)
                y = fmt % (float(toks[col["y"]]) / args.units_per_um)
                wf.write(f"{x}\t{y}\t{gene}\t{cnt}\t{toks[col['cell']]}\n")
                n_out += 1

    print(f"Converted {args.in_gem}: {n_in} rows read, {n_zero} dropped as zero-count, "
          f"{n_out} written, {len(seen_features)} distinct features -> {args.out}")

    # The cell counts are projected onto a model trained on the bin1 features; a naming
    # mismatch between the two GEMs produces empty cells rather than an error, so check.
    if args.check_features:
        want = read_check_features(args.check_features, args.check_features_colname)
        if not want:
            sys.exit(f"ERROR: no features read from {args.check_features} (--check-features).")
        overlap = len(want & seen_features) / len(want)
        if overlap < args.min_feature_overlap:
            examples_pixel = ", ".join(sorted(want)[:5])
            examples_cell = ", ".join(sorted(seen_features)[:5])
            # Drop the output: it is complete but useless, and leaving it behind invites
            # a later step to consume it as if the check had passed.
            os.remove(args.out)
            sys.exit(f"ERROR: only {overlap:.1%} of the {len(want)} pixel-level features appear in "
                     f"{args.in_gem} (need {args.min_feature_overlap:.0%}). The two GEMs most likely "
                     f"use different feature naming, which would silently yield empty cells.\n"
                     f"       pixel features: {examples_pixel}\n"
                     f"       cell features:  {examples_cell}\n"
                     f"       Set --colname-feature (geneName or geneID) so both ingests use the same column.")
        print(f"Feature check: {overlap:.1%} of {len(want)} pixel-level features present.")


if __name__ == "__main__":
    convert_stereoseq_cellbin(sys.argv[1:])
