suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(RcppParallel)
})
# if (!requireNamespace("this.path", quietly = TRUE)) {
#     install.packages("this.path", repos = "https://cran.rstudio.com")
# }

#script_path <- normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
#script_dir <- dirname(script_path)
script_dir <- this.path::this.dir()
source(file.path(script_dir, "umap_utils.r"))

parser <- ArgumentParser(description = "UMAP embedding on topic matrices")
parser$add_argument("--input",               type = "character", required = TRUE,
                    help = "Path to input TSV with metadata columns followed by numeric topic columns")
parser$add_argument("--out-prefix",                 type = "character", required = TRUE,
                    help = "Output file prefix (no extension)")
parser$add_argument("--tsv-colname-meta", nargs = "+", default = c("x", "y", "topK", "topP"),
                    help = "One or more column names for metadata in input TSV file")
parser$add_argument("--n-neighbors", type = "integer", default = NULL,
                    help = "UMAP n_neighbors. By default, it applies 50 when datasets has > 50 hexagons, 8 when <=50 and >20, 5 when <=20 and >10, 3 when <10.")
parser$add_argument("--pca-dims",            type = "integer",  default = 0,
                    help = "Optional PCA preprocessing dims for UMAP; 0 disables PCA (default: 0)")
parser$add_argument("--metric",              type = "character", default = "cosine",
                    help = "UMAP metric (default: cosine)")
parser$add_argument("--threads",             type = "integer",  default = 8,
                    help = "Number of threads for UMAP/RcppParallel (default: 8)")
parser$add_argument("--seed",                type = "integer",  default = 42,
                    help = "Random seed for reproducibility (default: 42)")
parser$add_argument("--no-sqrt", dest = "sqrt_transform", action = "store_false", default = TRUE,
                    help = "Disable sqrt transform on topic columns before UMAP")
args <- parser$parse_args()

set.seed(args$seed)
if (args$threads < 1) {
  stop("--threads must be >= 1")
}
setThreadOptions(numThreads = args$threads)

if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input)
}

log_message("Analysis started")
log_message(sprintf("Reading Input: %s", args$input))

df_lda <- fread(args$input, sep = "\t", header = TRUE)
log_message(sprintf("Metadata columns: %s", paste(args$tsv_colname_meta, collapse = ", ")))

matrix_info <- prepare_topic_matrix(df_lda, args$tsv_colname_meta, sqrt_transform = args$sqrt_transform)
mat <- matrix_info$matrix
log_message(sprintf("Matrix dims (after exclusion): %d rows x %d cols", nrow(mat), ncol(mat)))

# define the n-neighbors parameter based on dataset size
# Automatically adjust for very small datasets (≤ 50 rows)
pick_nn <- function(n, requested_k) {
  # If the user provides a valid k, clamp it within [2, n-1]
  if (!is.null(requested_k) && is.finite(requested_k) && requested_k >= 1) {
    return(max(2L, min(as.integer(requested_k), n - 1L)))
  }

  # Heuristic defaults for auto-selection
  if (n <= 10) {
    k <- 3L
  } else if (n <= 20) {
    k <- 5L
  } else if (n <= 50) {
    k <- 8L
  } else {
    k <- 50L
  }

  max(2L, min(k, n - 1L))
}

n_points <- nrow(mat)
k_final <- pick_nn(n_points, args$n_neighbors)
req <- args$n_neighbors
req_valid <- !is.null(req) && is.finite(req) && req >= 1

if (!req_valid) {
  # User didn't supply a valid k (NULL, NA, or <1): auto-selected via heuristic
  log_message(sprintf(
    "Auto-selected n_neighbors=%d for n=%d (requested=%s)",
    k_final, n_points,
    ifelse(is.null(req), "NULL", as.character(req))
  ))
} else if (k_final != as.integer(req)) {
  # User supplied a value but it was clamped to [2, n-1]
  log_message(sprintf(
    "Adjusted n_neighbors=%d for n=%d (requested=%s)",
    k_final, n_points, as.character(req)
  ))
} else {
  # Using the requested value as-is
  log_message(sprintf("Using requested n_neighbors=%d for n=%d", k_final, n_points))
}

log_message("Running UMAP (cosine)...")
umap_xy <- run_umap_embedding(
  mat,
  n_neighbors = k_final,
  threads     = args$threads,
  pca_dims    = args$pca_dims,
  metric      = args$metric
)

umap_df <- as.data.table(umap_xy)
setnames(umap_df, c("UMAP1", "UMAP2"))

meta_keep <- intersect(names(df_lda), args$tsv_colname_meta)
if (length(meta_keep) == 0) {
  log_message("No metadata columns found; proceeding with UMAP only.")
  out_dt <- umap_df
} else {
  out_dt <- cbind(df_lda[, ..meta_keep], umap_df)
}

umap_out_tsv <- paste0(args$out_prefix, ".umap.tsv.gz")
fwrite(out_dt, umap_out_tsv, sep = "\t")
log_message(paste("Wrote:", umap_out_tsv))
log_message("Done.")
