suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

script_path <- normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
script_dir <- dirname(script_path)
cartloader_root <- dirname(script_dir)
default_cmap <- file.path(cartloader_root, "assets", "fixed_color_map_256.tsv")
source(file.path(script_dir, "umap_utils.r"))

parser <- ArgumentParser(description = "Draw one UMAP scatter plot for each factor")
parser$add_argument("--input",                type = "character",   required = TRUE,
                    help = "Path to TSV containing UMAP coordinates and metadata")
parser$add_argument("--mode",                 type="character",     default="binary",
                    help = paste0(
                        "How to color the hexagons/cells:",
                        "'binary' = highlight hexagons/cells that belong to the current factor;",
                        "'prob'   = color hexagons/cells by their probability for the factor;",
                        "'both'   = generate both.",
                        "Default: binary"))
parser$add_argument("--tsv-colname-factor",   type = "character",   default = "topK",
                    help = "Column name in input used to colour the scatter (default: topK)")
parser$add_argument("--tsv-colname-umap1",       type = "character",   default = "UMAP1",
                    help = "Column name for UMAP X coordinate (default: UMAP1)")
parser$add_argument("--tsv-colname-umap2",       type = "character",   default = "UMAP2",
                    help = "Column name for UMAP Y coordinate (default: UMAP2)")
parser$add_argument("--tsv-colname-prob",    type = "character",   default = "topP",
                    help = "Column name for probability (used when --mode is prob or both; default: topP)")
parser$add_argument("--out",                  type = "character",   required = TRUE,
                    help = "Output prefix")
parser$add_argument("--subtitle",             type = "character",   default = NULL,
                    help = "Optional identifier used as the subtitle for the generated plot")

cmapgrp <- parser$add_argument_group("Optional Color Map Parameters.")
cmapgrp$add_argument("--cmap",                 type = "character", default = NULL,
                    help = "TSV with colour map. If omitted, use a fixed color map from cartloader (used if --mode is binary or both)")
cmapgrp$add_argument("--cmap-colname-factor",  type = "character", default = "Name",
                    help = "Column name in colour map for factor labels (default: Name) (used if --mode is binary or both)")
cmapgrp$add_argument("--cmap-colname-hex",     type = "character", default = "Color_hex",
                    help = "Column name in colour map for hex values (default: Color_hex) (used if --mode is binary or both)")

cmapgrp$add_argument("--missing-colour",       type = "character", default = "#a1a1a1ff",
                    help = "Hex colour to use for missing factor levels (default: #a1a1a1ff)")
cmapgrp$add_argument("--prob-palette",        type = "character", default = "viridis",
                    help = "Viridis palette option for probability mode (e.g., viridis, magma, plasma, inferno, cividis, turbo). Default: viridis")

plotgrp <- parser$add_argument_group("Optional Plot Parameters")
plotgrp$add_argument("--dpi",           type = "integer",  default = 300,
                    help = "DPI for the generated plot (default: 300)")
plotgrp$add_argument("--base-dim",      type = "double", default = 2,
                    help = "Baseline inches per panel added before scaling for auto-computed dimension (default: 2)")
plotgrp$add_argument("--scale-factor",  type = "double", default = 0.1,
                    help = "Inches per data unit of span for auto-computed dimension (default: 0.1)")
plotgrp$add_argument("--min-dim",       type = "double", default = 2,
                    help = "Minimum size in inches per panel for auto-computed dimension (default: 2)")
plotgrp$add_argument("--max-dim",       type = "double", default = 15,
                    help = "Maximum size in inches per panel for auto-computed dimension (default: 15)")
plotgrp$add_argument("--plot-dim",      type = "double", default = NULL, 
                    help = "Manual defined single plot size in inches (used for BOTH width and height) because it uses a 1:1 aspect. If omitted, the size is auto-computed from data span.")

args <- parser$parse_args()

if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input)
}

log_message("Rendering per-factor UMAP plots")
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

if (args$mode != "prob"){
  cmap_dt <- load_colormap(
    cmap_path = args$cmap,
    name_col  = args$cmap_colname_factor,
    hex_col   = args$cmap_colname_hex,
    default_cmap_path = default_cmap
  )
}

if (args$mode != "prob"){
  plot_prep <- prepare_plot_data(
    plot_dt,
    factor_col    = args$tsv_colname_factor,
    cmap_dt       = cmap_dt,
    color_column  = "color_key",
    missing_colour= args$missing_colour
  )
}


if (!is.null(args$plot_dim)) {
  log_message(paste0("Use user-defined dimension : ", args$plot_dim))
  plot_dim <- args$plot_dim
}else{
  log_message("Computing dimension based on the input...")
  umap1_dat <- plot_dt[[args$tsv_colname_umap1]]
  umap2_dat <- plot_dt[[args$tsv_colname_umap2]]
  plot_dim  <- auto_plot_dimension(umap1_dat, umap2_dat, 
                                    base_dim = args$base_dim, 
                                    scale_factor = args$scale_factor, 
                                    min_dim = args$min_dim, 
                                    max_dim = args$max_dim)
}

## Build a single faceted plot per mode (facet = factor)

# Determine factor levels to iterate
factor_vals <- as.character(plot_dt[[args$tsv_colname_factor]])
factor_levels <- unique(factor_vals[!is.na(factor_vals) & nzchar(factor_vals)])
if (length(factor_levels) == 0) {
  stop("No factor levels found to plot.")
}
# Order levels: numeric if possible else alphabetical
suppressWarnings({ fac_nums <- as.numeric(factor_levels) })
if (all(!is.na(fac_nums))) {
  ordered_levels <- factor_levels[order(fac_nums)]
} else {
  ordered_levels <- sort(factor_levels)
}

# Grid shape to size output
grid_ncol <- ceiling(sqrt(length(ordered_levels)))
grid_nrow <- ceiling(length(ordered_levels) / grid_ncol)

# Minimal working copy
xcol <- args$tsv_colname_umap1
ycol <- args$tsv_colname_umap2
fcol <- args$tsv_colname_factor
pcol <- args$tsv_colname_prob

dt_core <- plot_dt[, .(x = get(xcol), y = get(ycol), topf = as.character(get(fcol)))]
if (pcol %in% names(plot_dt)) {
  dt_core[, prob := as.numeric(plot_dt[[pcol]])]
} else {
  dt_core[, prob := NA_real_]
}

# Expand rows per facet level efficiently
make_panels_dt <- function() {
  dt <- rbindlist(lapply(ordered_levels, function(lbl) {
    tmp <- copy(dt_core)
    tmp[, panel_factor := lbl]
    tmp[, is_highlight := topf == lbl]
    tmp[, prob_value := ifelse(is_highlight, prob, NA_real_)]
    tmp
  }), use.names = TRUE)
  # Ensure facet order follows numeric interpretation when applicable
  dt[, panel_factor := factor(panel_factor, levels = ordered_levels)]
  dt
}

if (args$mode %in% c("binary", "both")) {
  if (!exists("plot_prep")) {
    # In case mode was changed programmatically; ensure palette exists
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
  }
  palette <- plot_prep$color_values
  dt_panels <- make_panels_dt()
  # Assign per-point color (highlight vs background)
  dt_panels[, highlight_col := palette[as.character(panel_factor)]]
  dt_panels[is.na(highlight_col), highlight_col := args$missing_colour]
  dt_panels[, color_hex := ifelse(is_highlight, highlight_col, args$missing_colour)]

  p_bin <- ggplot2::ggplot(dt_panels, ggplot2::aes(x = x, y = y, color = color_hex)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.7, shape = 16) +
    ggplot2::scale_color_identity() +
    ggplot2::facet_wrap(~ panel_factor, ncol = grid_ncol) +
    ggplot2::coord_equal() +
    theme_00(base_size = 11) +
    ggplot2::labs(title = "UMAP", subtitle = args$subtitle, x = xcol, y = ycol) +
    ggplot2::theme(legend.position = "none", 
                    panel.grid.minor = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(size = 16),
                    plot.subtitle = ggplot2::element_text(size = 12))

  out_png <- paste0(args$out, ".umap.single.binary.png")
  ggplot2::ggsave(
    out_png,
    plot = p_bin,
    width = plot_dim * grid_ncol,
    height = plot_dim * grid_nrow,
    units = "in",
    dpi = args$dpi,
    bg = "white"
  )
  log_message(paste("Wrote:", out_png))
}

if (args$mode %in% c("prob", "both")) {
  if (!pcol %in% names(plot_dt)) {
    stop(sprintf("Required probability column '%s' is missing for prob mode.", pcol))
  }
  dt_panels <- make_panels_dt()

  p_prob <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = dt_panels,
      ggplot2::aes(x = x, y = y),
      color = args$missing_colour, size = 0.25, alpha = 0.15, shape = 16
    ) +
    ggplot2::geom_point(
      data = dt_panels[is_highlight == TRUE],
      ggplot2::aes(x = x, y = y, color = prob_value),
      size = 0.35, alpha = 0.9, shape = 16
    ) +
    ggplot2::scale_color_viridis_c(option = args$prob_palette, limits = c(0, 1), oob = scales::squish) +
    ggplot2::facet_wrap(~ panel_factor, ncol = grid_ncol) +
    ggplot2::coord_equal() +
    theme_00(base_size = 11) +
    ggplot2::labs(title = "UMAP", subtitle = args$subtitle, x = xcol, y = ycol, color = "Prob") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(size = 16),
                    plot.subtitle = ggplot2::element_text(size = 12))

  out_png <- paste0(args$out, ".umap.single.prob.png")
  ggplot2::ggsave(
    out_png,
    plot = p_prob,
    width = plot_dim * grid_ncol,
    height = plot_dim * grid_nrow,
    units = "in",
    dpi = args$dpi,
    bg = "white"
  )
  log_message(paste("Wrote:", out_png))
}

log_message("Done.")
