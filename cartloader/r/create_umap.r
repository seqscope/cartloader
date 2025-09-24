suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(RcppParallel)
})

script_path <- normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
script_dir <- dirname(script_path)
source(file.path(script_dir, "umap_utils.r"))

parser <- ArgumentParser(description = "UMAP embedding on topic matrices")
parser$add_argument("--input",               type = "character", required = TRUE,
                    help = "Path to input TSV with metadata columns followed by numeric topic columns")
parser$add_argument("--out",                 type = "character", required = TRUE,
                    help = "Output file prefix (no extension)")
parser$add_argument("--tsv-colname-meta", nargs = "+", default = c("x", "y", "topK", "topP"),
                    help = "One or more column names for metadata in input TSV file")
parser$add_argument("--n-neighbors",         type = "integer",  default = 50,
                    help = "UMAP n_neighbors (default: 50)")
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

log_message("Running UMAP (cosine)...")
umap_xy <- run_umap_embedding(
  mat,
  n_neighbors = args$n_neighbors,
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

umap_out_tsv <- paste0(args$out, ".umap.tsv.gz")
fwrite(out_dt, umap_out_tsv, sep = "\t")
log_message(paste("Wrote:", umap_out_tsv))
log_message("Done.")
