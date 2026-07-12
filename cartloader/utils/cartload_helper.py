"""Helpers shared by run_cartload2 and run_cartload2_multi.

Keeping these in one place avoids the two packagers drifting apart (e.g. the
colormap normalization and the UMAP -> PMTiles command must match so that
per-sample catalogs and the multi-sample catalog agree).
"""
import os
from cartloader.utils.color_helper import normalize_rgb


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
