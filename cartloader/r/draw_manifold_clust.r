suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
})

# if (!requireNamespace("this.path", quietly = TRUE)) {
#     install.packages("this.path", repos = "https://cran.rstudio.com")
# }

script_dir <- this.path::this.dir()
cartloader_root <- dirname(dirname(script_dir))
default_cmap <- file.path(cartloader_root, "assets", "fixed_color_map_256.tsv")
source(file.path(script_dir, "umap_utils.r"))

parser <- ArgumentParser(description = "Draw one UMAP scatter, coloured by factor")
parser$add_argument("--tsv-manifold",                type = "character",   required = TRUE,
                    help = "Path to TSV containing UMAP/TSNE coordinates")
parser$add_argument("--tsv-clust",                type = "character",   required = TRUE,
                    help = "Path to TSV containing Cluster assignments")
parser$add_argument("--tsv-colname-ids",   type = "character",   nargs = "+",  default = "cell_id",
                    help = "Column name in input used to colour the scatter (default: cell_id)")
parser$add_argument("--tsv-colname-clust",   type = "character",   default = "topK",
                    help = "Column name in input used to colour the scatter (default: topK)")
parser$add_argument("--tsv-colname-x",       type = "character",   default = "UMAP1",
                    help = "Column name for manifold X coordinate (default: UMAP1)")
parser$add_argument("--tsv-colname-y",       type = "character",   default = "UMAP2",
                    help = "Column name for manifold Y coordinate (default: UMAP2)")
parser$add_argument("--out",             type = "character",   required = TRUE,
                    help = "Output file name")
parser$add_argument("--out-tsv",             type = "character", 
                    help = "Output joined TSV file name")
parser$add_argument("--subtitle",               type = "character",   default = NULL,
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

if (!file.exists(args$tsv_manifold)) {
  stop("Input file not found: ", args$tsv_manifold)
}

if (!file.exists(args$tsv_clust)) {
  stop("Input file not found: ", args$tsv_clust)
}

log_message("Rendering UMAP plot")
log_message(sprintf("Reading Input: %s", args$tsv_manifold))
df_manifold <- fread(args$tsv_manifold, sep = "\t", header = TRUE)

log_message(sprintf("Reading Input: %s", args$tsv_clust))
df_clust <- fread(args$tsv_clust, sep = "\t", header = TRUE)

log_message("Joining two input files")
plot_dt <- left_join(df_manifold, df_clust, by = args$tsv_colname_ids)

if (!is.null(args$out_tsv)) {
  log_message(sprintf("Writing joined TSV to: %s", args$out_tsv))
  fwrite(plot_dt, file = args$out_tsv, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE, compress = "auto")
}

required_cols <- c(args$tsv_colname_x, args$tsv_colname_y, args$tsv_colname_clust)
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

## check if the name columns can be converted to integer
if ( all(plot_dt[[args$tsv_colname_clust]] %>% is.integer()) ) {
  cmap_dt[, "cmap_factor" := as.factor(as.integer(cmap_dt$cmap_factor))]
  plot_dt[, (args$tsv_colname_clust) := as.factor(as.integer(get(args$tsv_colname_clust)))]
} else {
  cmap_dt[, "cmap_factor" := as.character(cmap_dt$cmap_factor)]
  plot_dt[, (args$tsv_colname_clust) := as.character(get(args$tsv_colname_clust))]
}

## join cmap_dt with plot_dt
plot_prep <- left_join(plot_dt, cmap_dt, by = setNames("cmap_factor", args$tsv_colname_clust))
## set cmap_hex as missing_colour if cmap_hex is missing
plot_prep[ is.na(cmap_hex), cmap_hex := args$missing_colour]

if (!is.null(args$plot_dim)) {
  log_message(paste0("Use user-defined dimension : ", args$plot_dim))
  plot_dim <- args$plot_dim
} else{
  log_message("Computing dimension based on the input...")
  x_dat <- plot_prep[[args$tsv_colname_x]]
  y_dat <- plot_prep[[args$tsv_colname_y]]
#   plot_dim  <- auto_plot_dimension(x_dat, y_dat, 
#                                     base_dim = args$base_dim, 
#                                     scale_factor = args$scale_factor, 
#                                     min_dim = args$min_dim, 
#                                     max_dim = args$max_dim)
  plot_dim <- 8
}

log_message(paste0("Plot dimension: ", plot_dim))

n_clust <- length(unique(plot_prep$cmap_hex))
max_clust_name_length <- max(nchar(as.character(unique(plot_prep$cmap_hex)))) + 3
legend_nrow <- ceiling((max_clust_name_length * n_clust) / 150)
n_points <- nrow(plot_prep)
pt_size <- 1 / log(n_points) / log(n_points) * 50

#print(head(plot_prep))

col_df <- plot_prep %>% distinct(!!sym(args$tsv_colname_clust), cmap_hex)
color_map <- setNames(col_df$cmap_hex, col_df[[args$tsv_colname_clust]])

p <- ggplot(plot_prep, aes(x = !!sym(args$tsv_colname_x), y = !!sym(args$tsv_colname_y))) +
    geom_point(aes(color = !!sym(args$tsv_colname_clust)), size = pt_size, alpha = 0.35, shape = 16) +
    labs(title = element_blank(), subtitle = args$subtitle, x = args$tsv_colname_x, y = args$tsv_colname_y) +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.key.width = grid::unit(0.6, "lines"),
        legend.key.height = grid::unit(0.5, "lines")
    ) +
    scale_color_manual(values = color_map, na.value = "grey70", drop = TRUE) +
    guides(color = guide_legend(nrow = legend_nrow, byrow = TRUE, override.aes = list(size = 3, alpha = 1)))

## if output file name ends with png, save as png,
if (grepl(".pdf$", args$out)) {
    umap_pdf = args$out
    ggsave( umap_pdf, 
        plot = p, 
        width = plot_dim, 
        height = plot_dim, 
        units = "in",
        bg = "white",
        limitsize = FALSE)
    log_message(paste("Wrote:", umap_pdf))
} else {
    if (grepl(".png$", args$out)) {
        umap_png = args$out
    } else {
        log_message("Output file does not end with .png or .pdf. Defaulting to .png")
        umap_png = paste0(args$out, ".png")
    }
    ggsave( umap_png, 
        plot = p, 
        width = plot_dim, 
        height = plot_dim, 
        units = "in",
        dpi = args$dpi, 
        bg = "white",
        limitsize = FALSE)
    log_message(paste("Wrote:", umap_png))
}
log_message("Done.")
