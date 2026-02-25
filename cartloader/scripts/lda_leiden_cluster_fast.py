import sys, os, argparse, inspect, gzip
import logging
import numpy as np
import igraph as ig
import scipy.sparse as sp
import leidenalg

# Faster IO
import polars as pl

# Faster cosine prep
from sklearn.preprocessing import normalize

# Approximate kNN
from pynndescent import NNDescent

from cartloader.utils.utils import create_custom_logger


def parse_arguments(_args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        prog=f"cartloader {inspect.getframeinfo(inspect.currentframe()).function}",
        description="Perform Leiden Clustering on LDA fit results",
    )

    run_params = parser.add_argument_group("Run Options")
    run_params.add_argument("--log", action="store_true", default=False, help="Write log to file")
    run_params.add_argument("--log-suffix", type=str, default=".log", help="Suffix for log file")

    io_params = parser.add_argument_group("Input/Output")
    io_params.add_argument("--tsv", type=str, required=True, help="Input TSV (optionally .gz)")
    io_params.add_argument("--out", type=str, required=True, help="Output CSV/TSV (optionally .gz)")
    io_params.add_argument("--offset-data", type=int, default=3, help="1-based column offset where numeric data starts")
    io_params.add_argument("--colname-cluster", type=str, default="cluster", help="Cluster column name")
    io_params.add_argument("--n-neighbors", type=int, default=50, help="k for kNN graph")
    io_params.add_argument("--random-seed", type=int, default=42, help="Random seed")
    io_params.add_argument("--key-ids", type=str, nargs="+", default=["cell_id"], help="Key ID columns")
    io_params.add_argument(
        "--drop-columns",
        type=str,
        nargs="*",
        default=["topK", "topP"],
        help='Columns to drop from numeric data (default: ["topK", "topP"])',
    )
    io_params.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution")
    io_params.add_argument("--use-weighting", action="store_true", default=False, help="Use weighted edges")
    io_params.add_argument("--sep", type=str, default="\t", help="Output separator")
    io_params.add_argument("--n-jobs", type=int, default=8, help="Number of threads for ANN")

    # Optional speed/quality knob (non-breaking to add)
    io_params.add_argument(
        "--mutual-knn",
        action="store_true",
        default=False,
        help="Keep only mutual kNN edges (often smaller graph, faster Leiden)",
    )

    if len(_args) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(_args)

def _build_knn_adjacency_pynndescent(
    X: np.ndarray,
    k: int,
    n_jobs: int,
    seed: int,
    weighted: bool,
    mutual_knn: bool,
):
    """
    Build a sparse kNN adjacency using approximate NNDescent.
    Returns CSR adjacency (n x n), symmetrized, with optional mutual-kNN pruning.
    """
    n = X.shape[0]
    # NNDescent neighbor_graph includes self in most settings; request k+1 and drop later.
    index = NNDescent(
        X,
        n_neighbors=k + 1,
        metric="cosine",
        random_state=seed,
        n_jobs=n_jobs,
        low_memory=True,
    )

    nbrs, dists = index.neighbor_graph  # shapes: (n, k+1)
    nbrs = nbrs[:, 1:]   # drop self
    dists = dists[:, 1:] # drop self

    # Build CSR efficiently without constructing an explicit "rows" array.
    # Each row has exactly k entries.
    indices = nbrs.reshape(-1).astype(np.int32, copy=False)

    if weighted:
        data = (1.0 - dists.reshape(-1)).astype(np.float32, copy=False)  # similarity-like
    else:
        data = np.ones(indices.shape[0], dtype=np.float32)

    indptr = (np.arange(0, n * k + 1, k)).astype(np.int64, copy=False)

    A = sp.csr_matrix((data, indices, indptr), shape=(n, n))

    # Remove any accidental self-loops (should be gone already, but safe)
    A.setdiag(0)
    A.eliminate_zeros()

    # Symmetrize to undirected
    A = A.maximum(A.T)

    # Optional: keep only mutual edges (often sparser & faster + sometimes cleaner clusters)
    if mutual_knn:
        # Keeps entries present in both directions (after sym this is still meaningful due to original directed nature)
        # A.multiply(A.T) works when A is not perfectly symmetric; after maximum it *is* symmetric,
        # but mutual_knn still helps if you move it before sym. We can approximate mutual by rebuilding from directed:
        # Here, we apply a conservative prune by requiring edge weight > 0 in both directions in original directed form.
        # Best approach is to mutualize before sym, but we didn't keep directed A; so we approximate by pruning weak links:
        # For unweighted, this does little; for weighted, you may prefer a threshold.
        # We'll do a true mutual by reconstructing directed first:
        # (A_dir is implicit in CSR construction; we can rebuild it quickly.)
        A_dir = sp.csr_matrix((data, indices, indptr), shape=(n, n))
        A_dir.setdiag(0)
        A_dir.eliminate_zeros()
        A = A_dir.multiply(A_dir.T)
        A.eliminate_zeros()
        # After mutual, make symmetric (it already is, but keep consistent)
        A = A.maximum(A.T)

    return A


def lda_leiden_cluster_fast(_args):
    """
    Perform Leiden clustering on LDA fit results
    """
    args = parse_arguments(_args)

    outprefix = os.path.splitext(args.out)[0]
    log_path = outprefix + "_lda_leiden_cluster" + args.log_suffix if args.log else None
    logger = create_custom_logger(__name__, log_path)

    logger.info("Analysis Started")
    logger.info(f"Reading results file {args.tsv} (Polars)...")

    # Polars can read .gz directly
    df = pl.read_csv(args.tsv, separator="\t", schema_overrides={key: pl.String for key in args.key_ids})

    n_rows, n_cols = df.shape
    logger.info(f"Loaded table: {n_rows} rows x {n_cols} cols")

    # Output (IDs only)
    missing_keys = [c for c in args.key_ids if c not in df.columns]
    if missing_keys:
        raise ValueError(f"Missing key-id columns in input: {missing_keys}")

    df_out = df.select(args.key_ids)

    # Feature columns: start at offset-data (1-based), drop any requested drop-columns
    start_idx = max(args.offset_data - 1, 0)
    if start_idx >= len(df.columns):
        raise ValueError(f"--offset-data={args.offset_data} is beyond number of columns ({len(df.columns)})")

    candidate_cols = df.columns[start_idx:]
    drop_set = set(args.drop_columns or []) | set(args.key_ids or [])
    feature_cols = [c for c in candidate_cols if c not in drop_set]

    if len(feature_cols) == 0:
        raise ValueError("No feature columns selected. Check --offset-data, --drop-columns, and --key-ids.")

    logger.info(f"Using {len(feature_cols)} feature columns (from column index {start_idx+1}).")
    logger.info(f"Dropping columns (if present): {args.drop_columns}")

    # Convert to float32 numpy and apply sqrt in-place
    X = df.select(feature_cols).to_numpy()
    X = np.asarray(X, dtype=np.float32, order="C")
    #np.sqrt(X + 1e-10, out=X)

    # if np.isinf(X.data if hasattr(X, "data") else X).any():
    #     print("Error: X contains Infinity values!")

    # # Check for NaNs correctly
    if np.isnan(X).any():
        print("Error: X contains actual NaNs!")
        nan_indices = np.argwhere(np.isnan(X))
        if len(nan_indices) > 0:
            row, col = nan_indices[0]
            print(f"FOUND NaN at Row: {row}, Col: {col}")
        else:
            print("No NaNs found.")

    # if X.min() < 0:
    #     raise ValueError("Negative values found in normalized data; cosine distance may behave unexpectedly.")

    # Normalize for cosine distance to behave well and often speed ANN
    normalize(X, norm="l2", axis=1, copy=False)

    logger.info(f"Constructing approximate kNN graph with PyNNDescent (k={args.n_neighbors}, n_jobs={args.n_jobs})...")
    A = _build_knn_adjacency_pynndescent(
        X,
        k=args.n_neighbors,
        n_jobs=args.n_jobs,
        seed=args.random_seed,
        weighted=args.use_weighting,
        mutual_knn=args.mutual_knn,
    )

    # Free large matrix ASAP
    del X
    try:
        import gc
        gc.collect()
    except Exception:
        pass

    # Convert adjacency to igraph without duplicates: keep upper triangle only
    logger.info("Building igraph (unique undirected edges)...")
    Au = sp.triu(A, k=1).tocoo()
    n = A.shape[0]

    g = ig.Graph(n=n, directed=False)
    # add_edges accepts iterables; avoid building a huge Python list if possible
    g.add_edges(zip(Au.row.tolist(), Au.col.tolist()))

    if args.use_weighting:
        g.es["weight"] = Au.data.astype(np.float32, copy=False)

    logger.info(f"Graph built: {g.vcount()} vertices, {g.ecount()} edges")

    logger.info("Running Leiden algorithm...")
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight" if args.use_weighting else None,
        resolution_parameter=args.resolution,
        seed=args.random_seed,
    )

    clusters = np.asarray(partition.membership, dtype=np.int32)
    df_out = df_out.with_columns(pl.Series(name=args.colname_cluster, values=clusters))

    logger.info(f"Writing output to {args.out}...")
    # Polars auto-compresses if the filename ends with .gz
    if args.out.endswith(".gz"):
        logger.info("Output will be gzip-compressed.")
        with gzip.open(args.out, "wb") as f_out:
            df_out.write_csv(f_out, separator=args.sep)
    else:
        df_out.write_csv(args.out, separator=args.sep)

    # if args.out.endswith(".gz"):
    #     import gzip
    #     with gzip.open(args.out, "wt") as f:
    #         df_out.write_csv(f, separator=args.sep)
    # else:
    #     df_out.write_csv(args.out, separator=args.sep)

    logger.info("Analysis completed.")


if __name__ == "__main__":
    # Base file name without extension
    script_name = os.path.splitext(os.path.basename(__file__))[0]

    # Dynamically get the function based on the script name
    func = getattr(sys.modules[__name__], script_name)

    # Call the function with command line arguments
    func(sys.argv[1:])