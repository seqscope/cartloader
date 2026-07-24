"""Shared feature include/exclude semantics.

One definition of "is this feature kept?", used by every Python step that filters
features by name (feature_select, reformat_cosmx, convert_stereoseq_cellbin). Regexes
are applied with ``re.search`` (substring semantics) to match the two other engines a
pattern may reach: sge_convert's pandas ``str.contains`` and spatula's
``std::regex_search``. Note that punkst matches with ``std::regex_match`` (whole
string) — never hand a pattern straight to punkst if it was written for these.
"""

import gzip, re, sys


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_feature_names(path):
    """Read the first column of a feature file into an ordered, de-duplicated list.

    Accepts both a bare list of names and a counts TSV (e.g. multi.features.tsv);
    '#'-prefixed and blank lines are skipped.
    """
    names, seen = [], set()
    with _open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            name = line.split("\t")[0].strip()
            if name and name not in seen:
                seen.add(name)
                names.append(name)
    return names


def read_feature_rows(path):
    """Read (name, count-or-None) pairs from a feature file, preserving file order."""
    rows, seen = [], set()
    with _open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            toks = line.split("\t")
            name = toks[0].strip()
            if not name or name in seen:
                continue
            seen.add(name)
            rows.append((name, toks[1].strip() if len(toks) > 1 else None))
    return rows


def compile_regex(pattern, flag_name):
    if pattern is None or pattern == "":
        return None
    try:
        return re.compile(pattern)
    except re.error as e:
        sys.exit(f"ERROR: invalid regex for {flag_name} ({e}): {pattern}")


class FeatureFilter:
    """Keeps a feature when it passes every filter that was supplied.

    An include filter narrows (the feature must be in the list / match the regex); an
    exclude filter removes. They combine, and an exclusion wins over an inclusion.
    """

    def __init__(self, include_list=None, exclude_list=None, include_regex=None, exclude_regex=None):
        self.keep_set = set(read_feature_names(include_list)) if include_list else None
        self.drop_set = set(read_feature_names(exclude_list)) if exclude_list else set()
        self.include_re = compile_regex(include_regex, "--include-feature-regex")
        self.exclude_re = compile_regex(exclude_regex, "--exclude-feature-regex")
        self.active = bool(self.keep_set is not None or self.drop_set
                           or self.include_re or self.exclude_re)
        self.n_dropped = 0

    def keep(self, name):
        ok = not (
            (self.keep_set is not None and name not in self.keep_set)
            or name in self.drop_set
            or (self.include_re is not None and not self.include_re.search(name))
            or (self.exclude_re is not None and self.exclude_re.search(name))
        )
        if not ok:
            self.n_dropped += 1
        return ok

    @classmethod
    def from_args(cls, args):
        """Build from an argparse namespace carrying the standard flag names."""
        return cls(include_list=getattr(args, "include_feature_list", None),
                   exclude_list=getattr(args, "exclude_feature_list", None),
                   include_regex=getattr(args, "include_feature_regex", None),
                   exclude_regex=getattr(args, "exclude_feature_regex", None))

    @staticmethod
    def add_arguments(group, what):
        """Register the four standard feature-filter flags on an argparse group."""
        group.add_argument('--include-feature-list', type=str, default=None,
                           help=f'Path to a file listing the feature names (one per line) to include in {what}')
        group.add_argument('--exclude-feature-list', type=str, default=None,
                           help=f'Path to a file listing the feature names (one per line) to exclude from {what}')
        group.add_argument('--include-feature-regex', type=str, default=None,
                           help=f'Regex of feature names to include in {what}')
        group.add_argument('--exclude-feature-regex', type=str, default=None,
                           help=f'Regex of feature names to exclude from {what}')
