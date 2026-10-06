import sys, os, argparse, inspect

from cartloader.utils.feature_filter import read_feature_names, read_feature_rows, compile_regex

def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Resolve feature include/exclude lists and regexes into a single feature list file.\n"
                    "Used to restrict a FICTURE2 analysis to a subset of features without altering the "
                    "transcript TSV or the feature list that packaging reads.")
    parser.add_argument('--in-features', type=str, default=None,
                        help='Universe of features: a TSV whose first column is the feature name and whose '
                             'optional second column is a count (e.g. multi.union_features.tsv). Lines starting '
                             'with "#" are skipped. Required unless --mode include with --include-list.')
    parser.add_argument('--out', required=True, type=str, help='Path to the output feature list')
    parser.add_argument('--mode', type=str, default="include", choices=["include", "exclude"],
                        help='What the output represents: "include" writes the features to keep (default), '
                             '"exclude" writes the features to drop (--exclude-list unioned with the '
                             '--in-features entries matching --exclude-regex).')
    parser.add_argument('--include-list', type=str, default=None, help='File listing feature names to include (one per line)')
    parser.add_argument('--exclude-list', type=str, default=None, help='File listing feature names to exclude (one per line)')
    parser.add_argument('--include-regex', type=str, default=None, help='Regex of feature names to include')
    parser.add_argument('--exclude-regex', type=str, default=None, help='Regex of feature names to exclude')
    parser.add_argument('--allow-empty', action='store_true', default=False,
                        help='Write an empty output instead of failing when no feature survives the filters')

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args(_args)

def feature_select(_args):
    args = parse_arguments(_args)

    for flag, path in (("--in-features", args.in_features), ("--include-list", args.include_list),
                       ("--exclude-list", args.exclude_list)):
        if path is not None and not os.path.exists(path):
            sys.exit(f"ERROR: file not found: {path} ({flag})")

    include_re = compile_regex(args.include_regex, "--include-regex")
    exclude_re = compile_regex(args.exclude_regex, "--exclude-regex")

    if args.mode == "exclude":
        if args.include_list or include_re:
            sys.exit("ERROR: --mode exclude writes the features to drop, so --include-list/--include-regex "
                     "cannot be combined with it.")
        excluded = []
        seen = set()
        for name in (read_feature_names(args.exclude_list) if args.exclude_list else []):
            if name not in seen:
                seen.add(name)
                excluded.append(name)
        if exclude_re is not None:
            if args.in_features is None:
                sys.exit("ERROR: --exclude-regex with --mode exclude needs --in-features to enumerate the "
                         "features the regex matches.")
            for name in read_feature_names(args.in_features):
                if name not in seen and exclude_re.search(name):
                    seen.add(name)
                    excluded.append(name)
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "wt") as wf:
            for name in excluded:
                wf.write(f"{name}\n")
        print(f"Wrote {len(excluded)} features to exclude: {args.out}")
        return

    # include mode: start from the universe, then apply every filter that was given.
    # Without --in-features the include list itself is the universe, which lets a plain
    # include list be filtered by a regex with no dataset-wide feature file at hand.
    if args.in_features is not None:
        rows = read_feature_rows(args.in_features)
        universe_src = args.in_features
    elif args.include_list is not None:
        rows = read_feature_rows(args.include_list)
        universe_src = args.include_list
    else:
        sys.exit("ERROR: --mode include needs --in-features (or --include-list to use as the universe).")

    keep = set(read_feature_names(args.include_list)) if args.include_list else None
    drop = set(read_feature_names(args.exclude_list)) if args.exclude_list else set()

    selected = []
    for name, count in rows:
        if keep is not None and name not in keep:
            continue
        if name in drop:
            continue
        if include_re is not None and not include_re.search(name):
            continue
        if exclude_re is not None and exclude_re.search(name):
            continue
        selected.append((name, count))

    if not selected and not args.allow_empty:
        sys.exit(f"ERROR: no feature survived the filters (universe: {universe_src}, "
                 f"{len(rows)} features). Check --include-list/--exclude-list feature names against that "
                 f"file, and --include-regex/--exclude-regex.")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    # Counts are carried over only when every selected feature has one: a file mixing
    # 2-column and 1-column rows makes punkst's feature reader warn and drop the sums.
    has_counts = all(c is not None for _, c in selected)
    with open(args.out, "wt") as wf:
        if has_counts:
            wf.write("#feature\ttotal_count\n")
        for name, count in selected:
            wf.write(f"{name}\t{count}\n" if has_counts else f"{name}\n")
    print(f"Selected {len(selected)} of {len(rows)} features: {args.out}")

if __name__ == "__main__":
    feature_select(sys.argv[1:])
