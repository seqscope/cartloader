import os, sys
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


def _pmtiles_format_flags(args):
    """punkst tile-op flags selecting the point tile format (args.tile_format, MVT by
    default). The feature property is named `gene`, as on the pmpoint/tippecanoe routes."""
    fmt = getattr(args, "tile_format", None) or "MVT"
    return ["--write-mlt-pmtiles" if fmt == "MLT" else "--write-mvt-pmtiles",
            "--feature-field-name", "gene"]


def _gene_bin_flags(args, feature_count_tsv, gene_bin_info):
    # "fixed" packs genes (sorted by count) into exactly --bin-count bins of about equal
    # molecule totals, the same assignment as spatula assign-feature2bin on the other
    # backends. punkst's default ("adaptive") treats --n-gene-bins only as an upper limit
    # and sizes bins to ~1M molecules each.
    parts = ["--feature-count-file", shlex.quote(feature_count_tsv),
             "--n-gene-bins", str(args.bin_count),
             "--gene-bin-mode", "fixed"]
    if gene_bin_info:
        # A precomputed gene->bin assignment (gene/count/bin rows, the spatula
        # assign-feature2bin format) takes precedence over the TSV-derived bins in punkst,
        # so every sample sharing it bins genes identically.
        parts.extend(["--gene-bin-info", shlex.quote(gene_bin_info)])
    return parts


def _ext_col_flags(layout):
    """punkst tile-op flags carrying the tiled transcripts' extra columns (e.g. cell_id)."""
    parts = []
    for flag, key in (("--ext-col-ints", "ext_ints"), ("--ext-col-floats", "ext_floats"),
                      ("--ext-col-strs", "ext_strs")):
        if layout[key]:
            parts.extend([flag] + [shlex.quote(x) for x in layout[key]])
    return parts


def make_direct_pmtiles_cmd(args, ficture2bin, out_prefix, tiled_prefix, feature_count_tsv, join_pixel_bins, join_pixel_ids, layout, gene_bin_info=None):
    """punkst tile-op: annotate the tiled molecules with the decoded pixel factors and write
    max-zoom point PMTiles (<out_prefix>_all.pmtiles and <out_prefix>_bin<N>.pmtiles, every
    archive carrying the factor and extra columns) plus <out_prefix>.pmtiles_index.tsv and
    <out_prefix>.bin_counts.json. Gene bins are always given: without them punkst keeps
    only the genes of the decoded model, dropping the molecules of other genes."""
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
    ] + _pmtiles_format_flags(args) + [
        "--pmtiles-zoom", str(args.point_max_zoom),
    ] + _gene_bin_flags(args, feature_count_tsv, gene_bin_info) + [
        "--out", shlex.quote(out_prefix),
        "--threads", str(args.threads),
    ]
    if len(join_pixel_bins) > 1:
        parts.extend(["--merge-emb"] + [shlex.quote(f"{x}.bin") for x in join_pixel_bins[1:]])
        parts.append("--merge-keep-all")
    if join_pixel_ids:
        parts.extend(["--emb-prefix"] + [shlex.quote(x) for x in join_pixel_ids])
    parts.extend(_ext_col_flags(layout))
    return " ".join(parts)


def make_dummy_factor_cmds(ficture2bin, dummy_prefix):
    """Commands writing a one-row factor source (<dummy_prefix>.bin/.index) far outside any
    dataset. punkst tile-op only writes point PMTiles through --annotate-pts, which needs
    some factor source; annotating against this one yields molecules with no factor
    columns (gene and count only)."""
    d = shlex.quote(os.path.dirname(dummy_prefix))
    tsv = shlex.quote(dummy_prefix + ".src.tsv")
    return [
        f"mkdir -p {d}",
        f"printf '#x\\ty\\tK1\\tP1\\n-1000000.00\\t-1000000.00\\t0\\t1.0000\\n' > {tsv}",
        " ".join([shlex.quote(ficture2bin), "pts2tiles", "--in-tsv", tsv,
                  "--tile-op-factor-tsv", "--binary-out", "--K 1", "--pixel-res 0.5",
                  "--tile-size 500", "--temp-dir", d, "--out-prefix", shlex.quote(dummy_prefix)]),
    ]


def make_compact_pmtiles_cmd(args, ficture2bin, out_prefix, tiled_prefix, feature_count_tsv, dummy_prefix, layout, gene_bin_info=None):
    """punkst tile-op against the dummy factor source: max-zoom point PMTiles of the tiled
    transcripts, carrying their columns (gene, count and extra columns such as cell_id) but
    no factor columns, with the same gene bins, index and bin-count JSON as
    make_direct_pmtiles_cmd given the same binning inputs."""
    parts = [
        shlex.quote(ficture2bin), "tile-op",
        "--in", shlex.quote(dummy_prefix),
        "--binary",
        "--annotate-pts", shlex.quote(tiled_prefix),
        "--icol-x", str(layout["idx_x"]),
        "--icol-y", str(layout["idx_y"]),
        "--icol-feature", str(layout["idx_feature"]),
        "--icol-count", str(layout["idx_count"]),
        "--anno-keep-all",
    ] + _pmtiles_format_flags(args) + [
        "--pmtiles-zoom", str(args.point_max_zoom),
    ] + _gene_bin_flags(args, feature_count_tsv, gene_bin_info) + [
        "--out", shlex.quote(out_prefix),
        "--threads", str(args.threads),
    ] + _ext_col_flags(layout)
    return " ".join(parts)


def make_index_filter_cmd(out_dir, index_tsv, drop_all=False):
    """Drop index rows whose PMTiles file was not written (punkst lists a gene bin that has
    no molecules in this sample but writes no archive for it), and optionally the `all`
    row, so neither the pyramid step nor the catalog refers to a missing file."""
    index_q = shlex.quote(index_tsv)
    out_dir_q = shlex.quote(out_dir)
    keep = '$1 != "all" && ' if drop_all else ""
    return (
        f"awk -F'\\t' -v d={out_dir_q} 'NR==1 {{ print; next }} "
        f"({keep}system(\"test -f \\\"\" d \"/\" $4 \"\\\"\") == 0) {{ print }}' "
        f"{index_q} > {index_q}.tmp && mv {index_q}.tmp {index_q}"
    )


def make_direct_pmtiles_pyramid_cmd(args, ficture2bin, out_dir, index_tsv, n_jobs=5):
    """Build the lower zoom levels of every point PMTiles listed in the index with
    `punkst build-pyramid`, n_jobs archives at a time, replacing each archive in place."""
    out_dir_q = shlex.quote(out_dir)
    index_q = shlex.quote(index_tsv)
    tmp_dir_q = shlex.quote(args.tmp_dir)
    punkst_q = shlex.quote(ficture2bin)
    xargs_opt = "-S 65536" if sys.platform == "darwin" else ""

    return (
        # Export python variables to bash environment variables to avoid nested quote issues
        f"export OUT_DIR={out_dir_q} TMP_DIR={tmp_dir_q} PUNKST={punkst_q}; "
        f"tail -n +2 {index_q} | cut -f4 | "
        # {{}} escapes the brackets in Python's f-string so xargs sees {}
        f"xargs -P {n_jobs} {xargs_opt} -I {{}} bash -c '"
        'rel="{}"; '
        'src="$OUT_DIR/$rel"; '
        'stem=$(basename "$src" .pmtiles); '
        # punkst keeps its temporary files next to --out and names the archive after it,
        # so build <tmp>/pyramid_<stem>/<stem>.pmtiles and move it over the input
        'tmpd="$TMP_DIR/pyramid_${stem}"; '
        'out="$tmpd/$stem.pmtiles"; '
        'mkdir -p "$tmpd" && '
        '"$PUNKST" build-pyramid --point-in "$src" --out "$out" '
        f"--min-zoom {args.point_min_zoom} "
        f"--max-tile-bytes {args.max_point_tile_bytes} "
        f"--max-tile-features {args.max_point_feature_counts} "
        f"--scale-factor-compression {args.pmpoint_compression_scale} "
        f"--threads {args.threads} && "
        # nothing is written when the input already reaches --min-zoom
        'if [ -f "$out" ]; then mv "$out" "$src"; fi && '
        'rmdir "$tmpd" 2>/dev/null || true'
        "'"
    )
