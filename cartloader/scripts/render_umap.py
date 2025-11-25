import argparse
import csv
import gzip
import json
import os
import sys
from typing import Iterable, List, Optional


def parse_arguments(_args):
    parser = argparse.ArgumentParser(prog="cartloader render_umap", description="Convert a UMAP TSV (optionally gzipped) into newline-delimited GeoJSON ready for tippecanoe.")
    parser.add_argument("--input", required=True, help="Path to the input TSV/TSV.GZ file")
    parser.add_argument("--out", required=True, help="Output NDJSON file")
    parser.add_argument("--delimiter", default="\t", help="Column delimiter for the input file (default: tab)")
    parser.add_argument("--colname-factor", default="topK", help="Column name containing factor IDs (falls back to --factor-column-candidates)")
    parser.add_argument("--colname-x", default="UMAP1", help="Column name containing the X coordinate (falls back to --x-column-candidates)")
    parser.add_argument("--colname-y", default="UMAP2", help="Column name containing the Y coordinate (falls back to --y-column-candidates)")
    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)
    
    args=parser.parse_args(_args)
    return args

# def _pick_column(headers: Iterable[str], preferred: Optional[str], fallbacks: Iterable[str]) -> Optional[str]:
#     header_set = {h: h for h in headers}
#     lower_map = {h.lower(): h for h in headers}

#     def _match(name: Optional[str]) -> Optional[str]:
#         if not name:
#             return None
#         if name in header_set:
#             return name
#         low = name.lower()
#         return lower_map.get(low)

#     candidate = _match(preferred)
#     if candidate:
#         return candidate
#     for item in fallbacks:
#         candidate = _match(item)
#         if candidate:
#             return candidate
#     return None

def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def render_umap(_args):
    args = parse_arguments(_args)

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    with _open_text(args.input) as fin:
        reader = csv.DictReader(fin, delimiter=args.delimiter)
        # headers = reader.fieldnames or []
        # x_col = _pick_column(headers, args.x_column, ["UMAP1", "umap1", "UMAP_1", "x", "X"] )
        # y_col = _pick_column(headers, args.y_column, ["UMAP2", "umap2", "UMAP_2", "y", "Y"])
        # factor_col = _pick_column(headers, args.factor_column, ["topK", "TopK", "factor", "Factor"])
        x_col =  args.colname_x
        y_col = args.colname_y
        factor_col = args.colname_factor

        missing = [ name for name, selected in (("X", x_col), ("Y", y_col), ("factor", factor_col)) if selected is None]
        if missing:
            raise SystemExit(f"Missing required columns in {args.input}: {', '.join(missing)}")

        total = 0
        with open(args.out, "w", encoding="utf-8") as fout:
            for row in reader:
                try:
                    x = float(row[x_col])  # type: ignore[arg-type]
                    y = float(row[y_col])  # type: ignore[arg-type]
                except (KeyError, ValueError, TypeError):
                    continue
                fid = row.get(factor_col)
                if not fid and fid != 0:
                    continue
                feature = {
                    "type": "Feature",
                    "geometry": {"type": "Point", "coordinates": [x, y]},
                    "properties": {"factor_id": str(fid)},
                }
                fout.write(json.dumps(feature, separators=(",", ":")))
                fout.write("\n")
                total += 1

        if total == 0:
            raise SystemExit("No valid rows were found while converting the UMAP TSV")

if __name__ == "__main__":
    # Get the base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])
