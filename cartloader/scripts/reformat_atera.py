import sys, os, argparse, gzip, inspect, datetime, shutil
import multiprocessing

import numpy as np

from cartloader.utils.feature_filter import FeatureFilter

UINT32_MAX = np.iinfo(np.uint32).max
# 10x cell ids print the uint32 prefix as 8 hex digits shifted to a-p ("0181d112" -> "abibnbbc")
HEX2CELLID = str.maketrans("0123456789abcdef", "abcdefghijklmnop")


def custom_log(message):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", file=sys.stderr)


def import_zarr():
    try:
        import zarr
    except ImportError:
        sys.exit("ERROR: reformat_atera reads the Atera *.zarr.zip files with zarr>=3, which is not "
                 "installed. Install it with: pip install 'zarr>=3'  (or pip install -e '.[atera]')")
    return zarr


def open_zip(path):
    """Open a *.zarr.zip in place (no unzipping); returns (store, root group)."""
    zarr = import_zarr()
    if not os.path.exists(path):
        sys.exit(f"ERROR: file not found: {path}")
    store = zarr.storage.ZipStore(path, mode="r")
    return store, zarr.open_group(store, mode="r")


def cell_id_strings(packed):
    """10x cell id strings for packed int64 ids (prefix + suffix << 32), e.g. 'abibnbbc-1'.
    -1 (no cell) becomes UNASSIGNED, the sentinel the cell stages already ignore."""
    out = []
    for p in packed.tolist():
        if p < 0:
            out.append("UNASSIGNED")
        else:
            out.append(f"{p & 0xFFFFFFFF:08x}".translate(HEX2CELLID) + f"-{p >> 32}")
    return out


def pack_cell_ids(raw):
    """(N, 2) uint32 [prefix, dataset suffix] -> int64 prefix + (suffix << 32); both uint32 max = -1."""
    packed = raw[:, 0].astype(np.int64) + (raw[:, 1].astype(np.int64) << 32)
    packed[(raw[:, 0] == UINT32_MAX) & (raw[:, 1] == UINT32_MAX)] = -1
    return packed


def parse_arguments(_args):
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Reformat a 10x Atera output bundle (transcripts/cells/cell_feature_matrix .zarr.zip, "
                    "read in place) into the generic transcript TSV + cell centroid (xy) + cell boundary "
                    "CSVs, plus one cluster CSV per clustering / cell annotation level in the bundle. "
                    "Coordinates stay in microns.")
    parser.add_argument('--in-dir', type=str, required=True, help='Atera output directory (containing transcripts.zarr.zip, cells.zarr.zip, cell_feature_matrix.zarr.zip)')
    parser.add_argument('--out', type=str, required=True, help='Output prefix for reformatted files')
    parser.add_argument('--threads', type=int, default=4, help='Worker processes converting transcript chunks in parallel (default: 4)')
    parser.add_argument('--chunk-rows', type=int, default=2000000, help='Transcripts per conversion chunk; bounds per-worker memory (default: 2000000)')
    parser.add_argument('--min-phred-score', type=float, default=20, help='Drop transcripts whose calibrated Q-score is below this (default: 20, as for Xenium)')
    parser.add_argument('--precision', type=int, default=2, help='Decimal places of the micron coordinates in the transcript TSV (default: 2)')
    parser.add_argument('--skip-transcripts', action='store_true', default=False, help='Write only the cell, boundary and cluster files')
    parser.add_argument('--out-suffix-tx', type=str, default='.transcripts.tsv.gz', help='Suffix for the output transcript TSV (default: .transcripts.tsv.gz)')
    parser.add_argument('--out-suffix-cells', type=str, default='.cells.csv.gz', help='Suffix for the output cell centroid CSV (default: .cells.csv.gz)')
    parser.add_argument('--out-suffix-boundaries', type=str, default='.cell_boundaries.csv.gz', help='Suffix for the output cell boundary CSV (default: .cell_boundaries.csv.gz)')
    parser.add_argument('--out-suffix-clusters', type=str, default='.clusters', help='Suffix stem for the cluster outputs: <out><stem>.<name>.csv and <out><stem>.<name>.alias.tsv (default: .clusters)')
    # Feature filtering happens as the transcript TSV is written, so a dropped feature is
    # absent from every downstream product of this ingest.
    ftr = parser.add_argument_group("Feature Filtering Parameters")
    FeatureFilter.add_arguments(ftr, "the output transcript TSV")
    return parser.parse_args(_args)


# ---------------------------------------------------------------------------
# Transcripts
# ---------------------------------------------------------------------------

# Per-worker state, set once by _init_worker: each process keeps its own zip handle open
# rather than re-reading the (large) zip directory for every chunk.
_W = {}


def _init_worker(tx_zip, gene_names, keep_gene, min_qv, precision):
    # One polars thread per worker: --threads already sets the parallelism, and polars'
    # default (one thread per core) in every worker oversubscribes the machine.
    os.environ["POLARS_MAX_THREADS"] = "1"
    import polars as pl
    _W["store"], _W["root"] = open_zip(tx_zip)
    _W["genes"] = pl.Series("gene", gene_names, dtype=pl.Utf8)
    _W["keep_gene"] = keep_gene
    _W["min_qv"] = min_qv
    _W["precision"] = precision


def _convert_chunk(task):
    """Convert rows [start, stop) of one spatial tile into a gzipped TSV part (no header)."""
    import polars as pl
    tile, start, stop, part = task
    t = _W["root"][f"grid/{tile}"]
    loc = np.asarray(t["location"][start:stop])
    qv = np.asarray(t["quality_score"][start:stop]).reshape(-1)
    cid = np.asarray(t["cell_id"][start:stop])
    nuc = np.asarray(t["overlaps_nucleus"][start:stop]).reshape(-1)
    gene_offset = np.asarray(t["gene_offset"][:])

    # Rows within a tile are sorted by gene and carry no gene column: gene_offset[g] is the
    # [start, end) row range of gene g, so row r belongs to the first gene whose range ends after r.
    gidx = np.searchsorted(gene_offset[:, 1], np.arange(start, stop), side="right")
    pass_qv = qv >= _W["min_qv"]
    keep = pass_qv & _W["keep_gene"][gidx]

    packed = pack_cell_ids(cid[keep])
    uniq, inv = np.unique(packed, return_inverse=True)
    df = pl.DataFrame({
        "X": loc[keep, 0],
        "Y": loc[keep, 1],
        "gene": _W["genes"].gather(gidx[keep]),
        "count": np.ones(int(keep.sum()), dtype=np.uint8),
        "cell_id": pl.Series(cell_id_strings(uniq), dtype=pl.Utf8).gather(inv.reshape(-1)),
        "Z": loc[keep, 2],
        "overlaps_nucleus": nuc[keep],
    })
    with gzip.open(part, "wb", compresslevel=3) as wf:
        wf.write(df.write_csv(separator="\t", include_header=False,
                              float_precision=_W["precision"]).encode())
    return stop - start, int(pass_qv.sum()), int(keep.sum())


def write_transcripts(args, gene_names, out_tsv):
    ftr_filter = FeatureFilter.from_args(args)
    keep_gene = np.array([ftr_filter.keep(g) for g in gene_names], dtype=bool)
    custom_log(f"  Keeping {keep_gene.sum()} of {len(gene_names)} features after feature filtering")

    tx_zip = os.path.join(args.in_dir, "transcripts.zarr.zip")
    store, root = open_zip(tx_zip)
    try:
        grid = root["grid"]
        tiles = [str(k) for k in grid.attrs["grid_keys"]]
        sizes = [grid[k]["location"].shape[0] for k in tiles]
    finally:
        store.close()
    custom_log(f"  {sum(sizes)} transcripts in {len(tiles)} tiles")

    part_dir = out_tsv + ".parts"
    os.makedirs(part_dir, exist_ok=True)
    tasks = []
    for tile, n in zip(tiles, sizes):
        for start in range(0, n, args.chunk_rows):
            tasks.append((tile, start, min(start + args.chunk_rows, n),
                          os.path.join(part_dir, f"part{len(tasks):06d}.tsv.gz")))

    with gzip.open(out_tsv, "wb") as wf:
        wf.write(b"X\tY\tgene\tcount\tcell_id\tZ\toverlaps_nucleus\n")
    # Concatenated gzip members form one valid gzip stream; imap keeps the task order, so
    # each part is appended (and deleted) as soon as it and every part before it are done.
    n_read = n_qv = n_kept = 0
    # spawn, not Linux's default fork: the parent has already run polars (cells/boundaries), and
    # a forked child inherits polars' thread-pool state without its threads, so it deadlocks.
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(args.threads, initializer=_init_worker,
              initargs=(tx_zip, gene_names, keep_gene, args.min_phred_score, args.precision)) as pool, \
         open(out_tsv, "ab") as out:
        for i, (nr, nq, nk) in enumerate(pool.imap(_convert_chunk, tasks)):
            part = tasks[i][3]
            with open(part, "rb") as pf:
                shutil.copyfileobj(pf, out)
            os.remove(part)
            n_read, n_qv, n_kept = n_read + nr, n_qv + nq, n_kept + nk
            if (i + 1) % 20 == 0 or i + 1 == len(tasks):
                custom_log(f"  Converted {i + 1}/{len(tasks)} chunks ({n_read} transcripts read)")
    os.rmdir(part_dir)
    custom_log(f"  Wrote {n_kept} transcripts to {out_tsv} ({n_read - n_qv} dropped by "
               f"Q-score < {args.min_phred_score:g}, {n_qv - n_kept} by feature filtering)")


# ---------------------------------------------------------------------------
# Cells, boundaries, clusters
# ---------------------------------------------------------------------------

def cluster_columns(attrs):
    """(output name, obs column) for each clustering and cell annotation level the bundle
    declares. The first clustering (graph-based Leiden) is named 'graphclust', Xenium's name
    for it and the one diffexp_genes/gene_expression_graphclust.csv uses; annotation levels
    keep their obs name (e.g. azimuth_fine)."""
    cols = []
    for i, c in enumerate(attrs.get("clusterings") or []):
        cols.append(("graphclust" if i == 0 else c["key_prefix"], c["key_prefix"]))
    for ann in attrs.get("cell_annotations") or []:
        for lvl in ann.get("cell_annotation_levels") or []:
            cols.append((lvl["level_prefix"], lvl["level_prefix"]))
    return cols


def write_clusters(obs, col, cell_ids, out_csv, out_alias):
    """Write <Barcode,Cluster> with dense 1-based integer ids (the cluster-file convention of
    the cell stages) and an index<TAB>alias label file keyed by the 0-based index they decode to.
    Only categories that hold at least one cell get an id, in category order, so no empty
    cluster shifts the indices."""
    codes = np.asarray(obs[f"{col}/codes"][:]).astype(np.int64)
    cats = [str(c) for c in obs[f"{col}/categories"][:]]
    used = np.unique(codes[codes >= 0])
    dense = np.zeros(len(cats), dtype=np.int64)
    dense[used] = np.arange(1, len(used) + 1)
    with gzip.open(out_csv, "wt") if out_csv.endswith(".gz") else open(out_csv, "w") as wf:
        wf.write("Barcode,Cluster\n")
        for cell, code in zip(cell_ids, codes.tolist()):
            if code >= 0:
                wf.write(f"{cell},{dense[code]}\n")
    with open(out_alias, "w") as wf:
        wf.write("index\talias\n")
        for i, code in enumerate(used.tolist()):
            wf.write(f"{i}\t{cats[code]}\n")
    custom_log(f"  {col}: {len(used)} clusters, {int((codes >= 0).sum())} cells -> {out_csv}")


def write_cells_and_clusters(args):
    """Cell centroids, boundaries and cluster files. Returns the feature names (transcript gene order)."""
    import polars as pl

    store, cells = open_zip(os.path.join(args.in_dir, "cells.zarr.zip"))
    try:
        packed = pack_cell_ids(np.asarray(cells["cell_id"][:]))
        cell_ids = np.array(cell_id_strings(packed), dtype=object)
        summary_cols = list(cells["cell_summary"].attrs["column_names"])
        summary = np.asarray(cells["cell_summary"][:])

        mstore, cfm = open_zip(os.path.join(args.in_dir, "cell_feature_matrix.zarr.zip"))
        try:
            obs = cfm["obs"]
            # obs shares cells.zarr's row order; checked, since cluster rows are matched by position
            obs_ids = np.asarray(obs["cell_id"][:]).astype(np.int64)
            if not np.array_equal(obs_ids, packed):
                sys.exit("ERROR: cell_feature_matrix obs rows are not in cells.zarr.zip cell order; "
                         "this Atera bundle layout is not supported yet.")
            tx_counts = np.asarray(obs["transcript_counts"][:])
            gene_names = [str(g) for g in cfm["var/feature_name"][:]]

            out_cells = f"{args.out}{args.out_suffix_cells}"
            custom_log(f"Writing {len(packed)} cell centroids to {out_cells}")
            with gzip.open(out_cells, "wb") as wf:
                wf.write(pl.DataFrame({
                    "cell_id": cell_ids.astype(str),
                    "x_centroid": summary[:, summary_cols.index("cell_centroid_x")],
                    "y_centroid": summary[:, summary_cols.index("cell_centroid_y")],
                    "transcript_counts": tx_counts,
                    "cell_area": summary[:, summary_cols.index("cell_area")],
                    "nucleus_area": summary[:, summary_cols.index("nucleus_area")],
                }).write_csv(float_precision=3).encode())

            custom_log("Writing cluster assignments")
            for name, col in cluster_columns(dict(cfm.attrs)):
                if col not in obs:
                    custom_log(f"  WARNING: obs has no '{col}' column; skipping")
                    continue
                stem = f"{args.out}{args.out_suffix_clusters}.{name}"
                write_clusters(obs, col, cell_ids, f"{stem}.csv", f"{stem}.alias.tsv")
        finally:
            mstore.close()

        # polygon_sets/1 = cell boundaries (0 = nuclei); vertices are fixed-width (x, y) rows
        # padded to the longest polygon, with num_vertices giving each row's valid pairs.
        ps = cells["polygon_sets/1"]
        n_poly = ps["num_vertices"].shape[0]
        out_bnd = f"{args.out}{args.out_suffix_boundaries}"
        custom_log(f"Writing {n_poly} cell boundaries to {out_bnd}")
        with gzip.open(out_bnd, "wb") as wf:
            wf.write(b"cell_id,vertex_x,vertex_y\n")
            for s in range(0, n_poly, 100000):
                e = min(s + 100000, n_poly)
                v = np.asarray(ps["vertices"][s:e])
                nv = np.asarray(ps["num_vertices"][s:e]).astype(np.int64)
                ci = np.asarray(ps["cell_index"][s:e]).astype(np.int64)
                max_v = v.shape[1] // 2
                xy = v.reshape(-1, max_v, 2)[np.arange(max_v)[None, :] < nv[:, None]]
                wf.write(pl.DataFrame({
                    "cell_id": cell_ids[np.repeat(ci, nv)].astype(str),
                    "vertex_x": xy[:, 0],
                    "vertex_y": xy[:, 1],
                }).write_csv(include_header=False, float_precision=3).encode())
    finally:
        store.close()
    return gene_names


def reformat_atera(_args):
    args = parse_arguments(_args)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    custom_log(f"Reading Atera bundle: {args.in_dir}")
    gene_names = write_cells_and_clusters(args)

    if args.skip_transcripts:
        custom_log("Skipping transcripts (--skip-transcripts)")
    else:
        out_tsv = f"{args.out}{args.out_suffix_tx}"
        custom_log(f"Writing transcripts to {out_tsv}")
        # Transcript rows index the cell-feature-matrix features; confirm the two agree.
        store, root = open_zip(os.path.join(args.in_dir, "transcripts.zarr.zip"))
        try:
            tx_genes = [str(g) for g in root.attrs["gene_names"]]
        finally:
            store.close()
        if tx_genes != gene_names:
            sys.exit("ERROR: transcripts.zarr.zip gene_names differ from cell_feature_matrix features.")
        write_transcripts(args, gene_names, out_tsv)
    custom_log("Done")


if __name__ == "__main__":
    reformat_atera(sys.argv[1:])
