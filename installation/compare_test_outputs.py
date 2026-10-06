#!/usr/bin/env python3
"""Check and compare the outputs of an end-to-end test run (used by docker_release.sh).

  compare_test_outputs.py check <out_dir>
      Every file referenced by every catalog.yaml under <out_dir> (the files upload_aws
      would upload) must exist and be non-empty. Exits 1 otherwise, or if there is no
      catalog.yaml at all.

  compare_test_outputs.py compare <new_out_dir> <baseline_out_dir> [--out report.md]
      Compares a run against a baseline run of the same test (e.g. with the previously
      released image) and writes a Markdown report: files missing from the new run, files
      only in the new run, line-count changes of text files, and large size changes. For
      review only; always exits 0, since outputs may legitimately differ between versions.

`check` needs cartloader (for upload_aws's catalog parsing), so run it inside the image.
"""
import argparse, gzip, os, sys

# Bookkeeping and scratch files that differ between runs by design.
SKIP_SUFFIXES = (".done", ".begin", ".mk", ".log")
SKIP_DIRS = {"tmp"}
TEXT_SUFFIXES = (".tsv", ".csv", ".txt", ".json", ".yaml", ".yml",
                 ".tsv.gz", ".csv.gz", ".txt.gz")
MAX_ROWS = 200

def check(out_dir):
    from cartloader.scripts.upload_aws import collect_files_from_yaml
    catalogs = sorted(os.path.join(d, f) for d, _, files in os.walk(out_dir)
                      for f in files if f == "catalog.yaml")
    if not catalogs:
        print(f"[FAIL] no catalog.yaml under {out_dir}")
        return 1
    n_bad = 0
    for catalog in catalogs:
        required, optional, basemaps = collect_files_from_yaml(catalog)
        files = sorted({v for v in required + optional + basemaps
                        if isinstance(v, str) and not v.startswith(("http://", "https://"))})
        base = os.path.dirname(catalog)
        bad = []
        for f in files:
            path = os.path.join(base, f)
            if not os.path.isfile(path):
                bad.append(f"missing: {f}")
            elif os.path.getsize(path) == 0:
                bad.append(f"empty:   {f}")
        status = "OK  " if not bad else "FAIL"
        print(f"[{status}] {catalog}: {len(files) - len(bad)}/{len(files)} referenced files present and non-empty")
        for b in bad:
            print(f"         {b}")
        n_bad += len(bad)
    return 1 if n_bad else 0

def list_files(root):
    out = {}
    for d, dirs, files in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
        for f in files:
            if f.endswith(SKIP_SUFFIXES):
                continue
            path = os.path.join(d, f)
            out[os.path.relpath(path, root)] = path
    return out

def count_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n

def pct(old, new):
    return f"{(new - old) / old * 100:+.1f}%" if old else "n/a"

def table(rows, header):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows[:MAX_ROWS]]
    if len(rows) > MAX_ROWS:
        lines.append(f"\n... and {len(rows) - MAX_ROWS} more")
    return "\n".join(lines)

def compare(new_dir, base_dir, out_path, size_tol):
    new, base = list_files(new_dir), list_files(base_dir)
    missing = sorted(set(base) - set(new))
    added = sorted(set(new) - set(base))
    common = sorted(set(new) & set(base))
    emptied, line_diff, size_diff = [], [], []
    for rel in common:
        s_new, s_base = os.path.getsize(new[rel]), os.path.getsize(base[rel])
        if s_new == 0 and s_base > 0:
            emptied.append(rel)
            continue
        if rel.endswith(TEXT_SUFFIXES):
            try:
                l_new, l_base = count_lines(new[rel]), count_lines(base[rel])
            except (OSError, EOFError) as e:
                line_diff.append((rel, "unreadable", str(e), ""))
                continue
            if l_new != l_base:
                line_diff.append((rel, l_base, l_new, pct(l_base, l_new)))
        if s_base and abs(s_new - s_base) / s_base > size_tol:
            size_diff.append((rel, s_base, s_new, pct(s_base, s_new)))

    md = [
        "# Test output comparison",
        "",
        f"- new: `{new_dir}`",
        f"- baseline: `{base_dir}`",
        "",
        f"{len(common)} files in both runs; **{len(missing)} missing from the new run**, "
        f"**{len(emptied)} empty only in the new run**, {len(added)} only in the new run, "
        f"{len(line_diff)} with different line counts, {len(size_diff)} with size changes "
        f"over {size_tol:.0%}.",
        "",
    ]
    if missing:
        md += ["## Missing from the new run", ""] + [f"- `{r}`" for r in missing[:MAX_ROWS]] + [""]
    if emptied:
        md += ["## Empty only in the new run", ""] + [f"- `{r}`" for r in emptied[:MAX_ROWS]] + [""]
    if added:
        md += ["## Only in the new run", ""] + [f"- `{r}`" for r in added[:MAX_ROWS]] + [""]
    if line_diff:
        md += ["## Line counts differ", "", table(line_diff, ["file", "baseline", "new", "change"]), ""]
    if size_diff:
        md += [f"## Size differs by more than {size_tol:.0%}", "",
               table(size_diff, ["file", "baseline bytes", "new bytes", "change"]), ""]
    text = "\n".join(md)
    if out_path:
        with open(out_path, "w") as fh:
            fh.write(text + "\n")
    print(text)
    return 0

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("check", help="Verify that every file referenced by the catalogs exists")
    p.add_argument("out_dir")
    p = sub.add_parser("compare", help="Compare a run against a baseline run (report only)")
    p.add_argument("new_out_dir")
    p.add_argument("baseline_out_dir")
    p.add_argument("--out", default=None, help="Write the Markdown report here as well")
    p.add_argument("--size-tolerance", type=float, default=0.2,
                   help="Report size changes larger than this fraction (default: 0.2)")
    args = parser.parse_args()
    if args.mode == "check":
        sys.exit(check(args.out_dir))
    sys.exit(compare(args.new_out_dir, args.baseline_out_dir, args.out, args.size_tolerance))

if __name__ == "__main__":
    main()
