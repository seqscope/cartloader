"""Helpers shared by run_cartload2 and run_cartload2_multi.

Keeping these in one place avoids the two packagers drifting apart (e.g. the
colormap normalization and the UMAP -> PMTiles command must match so that
per-sample catalogs and the multi-sample catalog agree).
"""
import os
import gzip
import json
from cartloader.utils.color_helper import normalize_rgb


def update_bin_counts_json(in_bin_json, in_features, out_bin_json, strip_comment_char="#"):
    """Rewrite a gene->bin assignment JSON so each gene's ``count`` reflects this
    sample's own feature totals, while keeping the shared ``bin`` assignment intact.

    When a shared _bin_counts.json (from spatula assign-feature2bin on the joint
    feature list) is reused across samples, its gene->bin routing must stay identical
    everywhere, but the per-gene ``count`` is the multi-sample aggregate. This reads
    the per-sample feature totals from ``in_features`` (column 0 = gene, column 1 =
    total count; comment/header rows are skipped) and replaces each record's ``count``
    with that sample's value. Genes absent from ``in_features`` get count 0.
    """
    counts = {}
    opener = gzip.open if in_features.endswith(".gz") else open
    with opener(in_features, "rt") as f:
        for line in f:
            if not line or line[0] == strip_comment_char:
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 2:
                continue
            try:
                counts[toks[0]] = int(toks[1])
            except ValueError:
                # non-numeric second column (e.g. a "gene<TAB>count" text header): skip
                continue

    with open(in_bin_json) as f:
        records = json.load(f)
    for rec in records:
        rec["count"] = counts.get(rec.get("gene"), 0)
    with open(out_bin_json, "w") as f:
        json.dump(records, f)


def copy_rgb_tsv(in_rgb, out_rgb, restart=False):
    """Normalize a FICTURE colormap (columns Name, R, G, B) into an rgb.tsv with
    columns Name, Color_index, R, G, B (RGB scaled to 0-1). Idempotent: only
    rewrites when missing, on --restart, or when the content would change."""
    def _get_content(path):
        with open(path, 'r') as f:
            hdrs = f.readline().rstrip().split("\t")
            col2idx = {hdr: i for i, hdr in enumerate(hdrs)}
            lines = ["\t".join(["Name", "Color_index", "R", "G", "B"])]
            for line in f:
                toks = line.rstrip().split("\t")
                if len(toks) != len(hdrs):
                    raise ValueError(f"Input RGB file {path} has inconsistent number of columns")
                rgb_r = float(toks[col2idx["R"]])
                rgb_g = float(toks[col2idx["G"]])
                rgb_b = float(toks[col2idx["B"]])
                rgb_r, rgb_g, rgb_b = normalize_rgb(rgb_r, rgb_g, rgb_b)
                name = toks[col2idx["Name"]]
                lines.append(f"{name}\t{name}\t{rgb_r:.5f}\t{rgb_g:.5f}\t{rgb_b:.5f}")
        return "\n".join(lines) + "\n"

    expected_content = _get_content(in_rgb)

    if restart or not os.path.exists(out_rgb):
        with open(out_rgb, 'w') as f:
            f.write(expected_content)
        return

    try:
        with open(out_rgb, 'r') as f:
            existing_content = f.read()
        if existing_content == expected_content:
            return  # up-to-date; no rewrite needed
    except Exception:
        pass

    with open(out_rgb, 'w') as f:
        f.write(expected_content)


def record_catalog_alias(catalog_path, factor_id, alias_filename, key="alias"):
    """Record ``key: alias_filename`` on the factor whose id is ``factor_id`` in a
    CartoScope catalog YAML, editing it in place.

    Handles both catalog shapes: ``assets.factors`` as a list of dicts (the
    per-sample catalog.yaml) and as a map keyed by hyphenated id (the
    multi-catalog.yaml). Raises if the factor id is not present, so a mistyped
    projection id fails the run rather than silently deploying an unreferenced file.
    """
    import yaml
    with open(catalog_path) as f:
        catalog = yaml.safe_load(f)
    factors = catalog.get("assets", {}).get("factors")
    if factors is None:
        factors = catalog.get("factors")  # legacy top-level map
    matched = False
    if isinstance(factors, dict):
        entry = factors.get(factor_id)
        if entry is not None:
            entry[key] = alias_filename
            matched = True
    elif isinstance(factors, list):
        for entry in factors:
            if entry.get("id") == factor_id:
                entry[key] = alias_filename
                matched = True
                break
    if not matched:
        raise ValueError(f"record_catalog_alias: factor '{factor_id}' not found in {catalog_path}")
    with open(catalog_path, "w") as f:
        yaml.dump(catalog, f, Dumper=yaml.SafeDumper, default_flow_style=False, sort_keys=False)


def render_umap_cmd(in_tsv, out_ndjson, colname_factor="topK", colname_x="UMAP1", colname_y="UMAP2"):
    """Command to convert a UMAP TSV into NDJSON points (via cartloader render_umap)."""
    return " ".join([
        "cartloader", "render_umap",
        f"--input {in_tsv}",
        f"--out {out_ndjson}",
        f"--colname-factor {colname_factor}",
        f"--colname-x {colname_x}",
        f"--colname-y {colname_y}",
    ])


def umap_tippecanoe_cmd(out_pmtiles, in_ndjson, tippecanoe, tmp_dir,
                        threads=4, min_zoom=0, max_zoom=18, preserve_thres=1024):
    """Command to build a UMAP PMTiles pyramid from NDJSON points (via tippecanoe)."""
    return " ".join([
        f"TIPPECANOE_MAX_THREADS={threads}",
        f"'{tippecanoe}'",
        f"-t {tmp_dir}",
        f"-o {out_pmtiles}",
        "-Z", str(min_zoom),
        "-z", str(max_zoom),
        "-l", "umap",
        "--force",
        "--drop-densest-as-needed",
        "--extend-zooms-if-still-dropping",
        "--no-duplication",
        f"--preserve-point-density-threshold={preserve_thres}",
        in_ndjson,
    ])
