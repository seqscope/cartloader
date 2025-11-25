suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

script_path <- normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
script_dir <- dirname(script_path)
cartloader_root <- dirname(script_dir)
default_cmap <- file.path(cartloader_root, "assets", "fixed_color_map_256.tsv")
source(file.path(script_dir, "umap_utils.r"))

parser <- ArgumentParser(description = "Draw one UMAP scatter, coloured by factor")
parser$add_argument("--input",                type = "character",   required = TRUE,
                    help = "Path to TSV containing UMAP coordinates and metadata")
parser$add_argument("--tsv-colname-factor",   type = "character",   default = "topK",
                    help = "Column name in input used to colour the scatter (default: topK)")
parser$add_argument("--tsv-colname-umap1",       type = "character",   default = "UMAP1",
                    help = "Column name for UMAP X coordinate (default: UMAP1)")
parser$add_argument("--tsv-colname-umap2",       type = "character",   default = "UMAP2",
                    help = "Column name for UMAP Y coordinate (default: UMAP2)")
parser$add_argument("--out-prefix",                  type = "character",   required = TRUE,
                    help = "Output prefix")
parser$add_argument("--subtitle",             type = "character",   default = NULL,
                    help = "Optional identifier used as the subtitle for the generated plot")

cmapgrp <- parser$add_argument_group("Optional Color Map Parameters")
cmapgrp$add_argument("--cmap",                 type = "character", default = NULL,
                    help = "TSV with colour map. If omitted, use a fixed color map from cartloader")
cmapgrp$add_argument("--cmap-colname-factor",  type = "character", default = "Name",
                    help = "Column name in colour map for factor labels (default: Name)")
cmapgrp$add_argument("--cmap-colname-hex",     type = "character", default = "Color_hex",
                    help = "Column name in colour map for hex values (default: Color_hex)")
cmapgrp$add_argument("--missing-colour",       type = "character", default = "#a1a1a1ff",
                    help = "Hex colour to use for missing factor levels (default: #a1a1a1ff)")

plotgrp <- parser$add_argument_group("Optional Plot Parameters")
plotgrp$add_argument("--dpi",           type = "integer",  default = 300,
                    help = "DPI for the generated plot (default: 300)")
plotgrp$add_argument("--base-dim",      type = "double", default = 6,
                 help = "Baseline inches added before scaling for auto-computed dimension (default: 6)")
plotgrp$add_argument("--scale-factor",  type = "double", default = 0.2,
                 help = "Inches per data unit of span for auto-computed dimension (default: 0.2)")
plotgrp$add_argument("--min-dim",       type = "double", default = 4,
                 help = "Minimum size in inches for auto-computed dimension (default: 4)")
plotgrp$add_argument("--max-dim",       type = "double", default = 15,
                 help = "Maximum size in inches for auto-computed dimension (default: 15)")
plotgrp$add_argument("--plot-dim",      type = "double", default = NULL, 
                  help = "Manual defined single plot size in inches (used for BOTH width and height) because it uses a 1:1 aspect. If omitted, the size is auto-computed from data span.")

args <- parser$parse_args()

if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input)
}

log_message("Rendering UMAP plot")
log_message(sprintf("Reading Input: %s", args$input))

plot_dt <- fread(args$input, sep = "\t", header = TRUE)
required_cols <- c(args$tsv_colname_umap1, args$tsv_colname_umap2, args$tsv_colname_factor)
missing_cols <- setdiff(required_cols, names(plot_dt))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Input is missing required columns: %s",
    paste(missing_cols, collapse = ", ")
  ))
}

cmap_dt <- load_colormap(
  cmap_path = args$cmap,
  name_col  = args$cmap_colname_factor,
  hex_col   = args$cmap_colname_hex,
  default_cmap_path = default_cmap
)

plot_prep <- prepare_plot_data(
  plot_dt,
  factor_col    = args$tsv_colname_factor,
  cmap_dt       = cmap_dt,
  color_column  = "color_key",
  missing_colour= args$missing_colour
)


if (!is.null(args$plot_dim)) {
  log_message(paste0("Use user-defined dimension : ", args$plot_dim))
  plot_dim <- args$plot_dim
}else{
  log_message("Computing dimension based on the input...")
  umap1_dat <- plot_prep$data[[args$tsv_colname_umap1]]
  umap2_dat <- plot_prep$data[[args$tsv_colname_umap2]]
  plot_dim  <- auto_plot_dimension(umap1_dat, umap2_dat, 
                                    base_dim = args$base_dim, 
                                    scale_factor = args$scale_factor, 
                                    min_dim = args$min_dim, 
                                    max_dim = args$max_dim)
}

p <- draw_umap_plot(
  plot_data    = plot_prep$data,
  color_values = plot_prep$color_values,
  color_column = plot_prep$color_column,
  subtitle     = args$subtitle,
  x_col        = args$tsv_colname_umap1,
  y_col        = args$tsv_colname_umap2
)

umap_png = paste0(args$out_prefix, ".umap.png")
ggsave( umap_png, 
        plot = p, 
        width = plot_dim, 
        height = plot_dim, 
        units = "in",
        dpi = args$dpi, 
        bg = "white",
        limitsize = FALSE)


log_message(paste("Wrote:", umap_png))
log_message("Done.")
