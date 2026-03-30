import os
import shlex


def infer_tiled_query_layout(tiled_prefix, colname_feature="gene", colname_count="count"):
    """
    Infer query-column indices and extra columns from a tiled transcript TSV.
    """
    tsv_path = tiled_prefix + ".tsv"
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"File not found: {tsv_path}")

    header = None
    sample_rows = []
    with open(tsv_path, "r") as f:
        header_buffer = []
        for raw in f:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("##"):
                header_buffer.append(line)
                continue
            if line.startswith("#"):
                header_buffer.append(line)
                continue
            break
        if header_buffer:
            header = header_buffer[-1].lstrip("#").split("\t")
    if header is None:
        raise ValueError(f"Cannot find a header line in tiled transcript TSV: {tsv_path}")

    with open(tsv_path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            sample_rows.append(line.split("\t"))
            if len(sample_rows) >= 100:
                break

    lower2idx = {name.lower(): idx for idx, name in enumerate(header)}

    def _find_idx(candidates, fallback):
        for name in candidates:
            if name.lower() in lower2idx:
                return lower2idx[name.lower()]
        return fallback

    idx_x = _find_idx(["x", "lon"], 0)
    idx_y = _find_idx(["y", "lat"], 1)
    idx_feature = _find_idx([colname_feature, "gene", "feature"], 2)
    idx_count = _find_idx([colname_count, "count", "ct", "gn"], 3)
    na_tokens = ["", ".", "NA", "NaN", "nan"]
    na_token_set = set(na_tokens)

    def _guess_extra_type(name, tokens):
        lname = name.lower()
        if "id" in lname:
            return "str"
        if lname in {"z", "z_location", "zloc"} or lname.startswith("z_"):
            return "float"

        saw_float = False
        saw_string = False
        saw_value = False
        for token in tokens:
            if token in na_token_set:
                continue
            saw_value = True
            try:
                int(token)
                continue
            except Exception:
                pass
            try:
                float(token)
                saw_float = True
                continue
            except Exception:
                saw_string = True
                break
        if saw_string:
            return "str"
        if saw_float:
            return "float"
        if saw_value:
            return "int"
        return "str"

    def _guess_null_token(tokens):
        na_counts = {}
        for token in tokens:
            if token in na_token_set:
                na_counts[token] = na_counts.get(token, 0) + 1
        if not na_counts:
            return None
        return max(na_tokens, key=lambda token: (na_counts.get(token, 0), -na_tokens.index(token)))

    ext_ints = []
    ext_floats = []
    ext_strs = []
    reserved = {idx_x, idx_y, idx_feature, idx_count}
    for idx, name in enumerate(header):
        if idx in reserved:
            continue
        tokens = [row[idx] for row in sample_rows if idx < len(row)]
        null_token = _guess_null_token(tokens)
        spec = f"{idx}:{name}:{null_token}" if null_token is not None else f"{idx}:{name}"
        typ = _guess_extra_type(name, tokens)
        if typ == "int":
            ext_ints.append(spec)
        elif typ == "float":
            ext_floats.append(spec)
        else:
            ext_strs.append(spec)

    return {
        "idx_x": idx_x,
        "idx_y": idx_y,
        "idx_feature": idx_feature,
        "idx_count": idx_count,
        "ext_ints": ext_ints,
        "ext_floats": ext_floats,
        "ext_strs": ext_strs,
    }


def make_direct_pmtiles_cmd(args, ficture2bin, out_prefix, tiled_prefix, feature_count_tsv, join_pixel_bins, join_pixel_ids, layout):
    parts = [
        shlex.quote(ficture2bin), "tile-op",
        "--in", shlex.quote(join_pixel_bins[0]),
        "--binary",
        "--annotate-pts", shlex.quote(tiled_prefix),
        "--icol-x", str(layout["idx_x"]),
        "--icol-y", str(layout["idx_y"]),
        "--icol-feature", str(layout["idx_feature"]),
        "--icol-count", str(layout["idx_count"]),
        "--anno-keep-all",
        "--write-mlt-pmtiles",
        "--pmtiles-zoom", str(args.point_max_zoom),
        "--feature-count-file", shlex.quote(feature_count_tsv),
        "--n-gene-bins", str(args.bin_count),
        "--out", shlex.quote(out_prefix),
        "--threads", str(args.threads),
    ]
    if len(join_pixel_bins) > 1:
        parts.extend(["--merge-emb"] + [shlex.quote(f"{x}.bin") for x in join_pixel_bins[1:]])
        parts.append("--merge-keep-all")
    if join_pixel_ids:
        parts.extend(["--emb-prefix"] + [shlex.quote(x) for x in join_pixel_ids])
    if layout["ext_ints"]:
        parts.extend(["--ext-col-ints"] + [shlex.quote(x) for x in layout["ext_ints"]])
    if layout["ext_floats"]:
        parts.extend(["--ext-col-floats"] + [shlex.quote(x) for x in layout["ext_floats"]])
    if layout["ext_strs"]:
        parts.extend(["--ext-col-strs"] + [shlex.quote(x) for x in layout["ext_strs"]])
    return " ".join(parts)


def make_direct_pmtiles_pyramid_cmd(args, out_dir, index_tsv):
    out_dir_q = shlex.quote(out_dir)
    index_q = shlex.quote(index_tsv)
    tmp_dir_q = shlex.quote(args.tmp_dir)
    pmpoint_q = shlex.quote(args.pmpoint)
    cleanup = "" if args.keep_intermediate_files else 'rm -rf "$tmpd"; '
    return (
        f"tail -n +2 {index_q} | cut -f4 | while read rel; do "
        f'src={out_dir_q}/"$rel"; '
        'stem=$(basename "$src" .pmtiles); '
        f'tmpd={tmp_dir_q}/"${{stem}}"; '
        'out="$src.tmp"; '
        'mkdir -p "$tmpd"; '
        f"{pmpoint_q} build-pyramid-pmtiles "
        f"--scale-factor-compression {args.pmpoint_compression_scale} "
        f'--tmp-dir "$tmpd" --in "$src" --out "$out" '
        f"--min-zoom {args.point_min_zoom} "
        f"--max-tile-bytes {args.max_point_tile_bytes} "
        f"--max-tile-features {args.max_point_feature_counts} "
        f"--threads {args.threads}; "
        'mv "$out" "$src"; '
        f"{cleanup}"
        "done"
    )
