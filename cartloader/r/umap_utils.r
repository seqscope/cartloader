suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(uwot)
})

log_message <- function(msg) {
  cat(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), "\n", sep = "")
}

prepare_topic_matrix <- function(df, meta_cols, sqrt_transform = TRUE, zero_fill = TRUE) {
  topic_cols <- setdiff(names(df), meta_cols)
  if (length(topic_cols) == 0) {
    stop("No topic columns detected after excluding metadata columns.")
  }

  topic_dt <- copy(df[, ..topic_cols])
  original_na_counts <- vapply(topic_dt, function(col) sum(is.na(col)), integer(1))
  topic_dt[, (topic_cols) := lapply(.SD, function(col) type.convert(col, as.is = TRUE)), .SDcols = topic_cols]

  non_numeric_cols <- topic_cols[!vapply(topic_dt, is.numeric, logical(1))]
  if (length(non_numeric_cols) > 0) {
    stop(sprintf("Topic columns must be numeric; unable to convert: %s", paste(non_numeric_cols, collapse = ", ")))
  }

  new_na_counts <- vapply(topic_dt, function(col) sum(is.na(col)), integer(1))
  na_diff <- new_na_counts - original_na_counts
  if (any(na_diff > 0)) {
    warn_cols <- names(na_diff)[na_diff > 0]
    log_message(sprintf(
      "Warning: Non-numeric values coerced to NA in columns: %s",
      paste(warn_cols, collapse = ", ")
    ))
  }

  mat <- as.matrix(topic_dt)

  if (anyNA(mat)) {
    n_na <- sum(is.na(mat))
    log_message(sprintf("Warning: %d NA values in numeric matrix; coercing NAs to 0.", n_na))
    if (isTRUE(zero_fill)) {
      mat[is.na(mat)] <- 0
    }
  }
  if (any(mat < 0, na.rm = TRUE)) {
    log_message("Warning: negative values found; clipping to 0 before sqrt.")
    mat[mat < 0] <- 0
  }

  if (isTRUE(sqrt_transform)) {
    mat <- sqrt(mat)
  }

  list(matrix = mat, topic_cols = topic_cols)
}

run_umap_embedding <- function(mat, n_neighbors = 50, threads = 8, pca_dims = 0, metric = "cosine") {
  pca_arg <- if (pca_dims > 0) pca_dims else NULL
  uwot::umap2(
    X            = mat,
    metric       = metric,
    n_components = 2,
    n_neighbors  = n_neighbors,
    n_threads    = threads,
    pca          = pca_arg
  )
}

load_colormap <- function(
  cmap_path = NULL,
  name_col = "Name",
  hex_col  = "Color_hex",
  default_cmap_path = NULL,
  strict_hex = TRUE
) {
  resolved_path <- cmap_path
  if ((is.null(resolved_path) || identical(resolved_path, "")) && !is.null(default_cmap_path)) {
    resolved_path <- default_cmap_path
    log_message(sprintf("Reading default color map: %s", resolved_path))
  } else if (!is.null(resolved_path) && !identical(resolved_path, "")) {
    log_message(sprintf("Reading user provided color map: %s", resolved_path))
  }

  if (is.null(resolved_path) || identical(resolved_path, "")) {
    stop("No colormap path provided and no default available.")
  }

  if (!file.exists(resolved_path)) {
    stop(sprintf("Colormap file not found: %s", resolved_path))
  }

  cmap_dt <- fread(resolved_path, sep = "\t", header = TRUE)
  required_cols <- c(name_col, hex_col)
  if (!all(required_cols %in% names(cmap_dt))) {
    stop(sprintf(
      "Colormap file must contain columns: %s",
      paste(required_cols, collapse = ", ")
    ))
  }

  cmap_dt <- copy(cmap_dt)
  cmap_dt[, cmap_factor := as.character(get(name_col))]
  cmap_dt[, cmap_hex := as.character(get(hex_col))]

  if (strict_hex) {
    invalid_hex <- cmap_dt[!is.na(cmap_hex) & !grepl("^#[0-9A-Fa-f]{6}$", cmap_hex), unique(cmap_hex)]
    if (length(invalid_hex) > 0) {
      stop(sprintf(
        "Invalid hex colour codes detected: %s",
        paste(invalid_hex, collapse = ", ")
      ))
    }
  }

  cmap_dt=cmap_dt[,c("cmap_factor", "cmap_hex")]

  setkey(cmap_dt, cmap_factor)
  cmap_dt
}

prepare_plot_data <- function(
  plot_dt,
  factor_col,
  cmap_dt,
  color_column = "color_key",
  missing_colour = "#BEBEBE"
) {
  if (!factor_col %in% names(plot_dt)) {
    stop(sprintf("Column '%s' is required in results to apply the colormap.", factor_col))
  }

  working_dt <- copy(plot_dt)
  # Detect true missing before normalization to control legend entry
  missing_mask <- is.na(working_dt[[factor_col]]) | working_dt[[factor_col]] == ""
  n_missing_vals <- sum(missing_mask, na.rm = TRUE)
  working_dt[, (factor_col) := as.character(get(factor_col))]
  working_dt[missing_mask, (factor_col) := "(missing)"]

  setkeyv(working_dt, factor_col)

  # joined_dt <- cmap_dt[working_dt, on = setNames(factor_col, "cmap_factor")]
  # joined_dt <- cmap_dt[working_dt, on = setNames("cmap_factor", factor_col)]
  joined_dt <- merge(
      working_dt,
      cmap_dt,
      by.x = factor_col,
      by.y = "cmap_factor",
      all.x = TRUE,
      sort = FALSE
    )

  # Assign fallback colour to any rows without a mapped hex
  joined_dt[is.na(cmap_hex), cmap_hex := missing_colour]

  # Base palette from cmap file
  base_palette <- stats::setNames(as.character(cmap_dt$cmap_hex), as.character(cmap_dt$cmap_factor))
  base_palette <- base_palette[!is.na(names(base_palette)) & names(base_palette) != "" & !is.na(base_palette)]

  # Determine levels actually used in the data
  used_levels <- as.character(sort(unique(joined_dt[[factor_col]])))

  # Add any missing used levels (no cmap mapping) with fallback colour
  missing_levels <- setdiff(used_levels, names(base_palette))
  if (length(missing_levels) > 0) {
    log_message(sprintf(
      "Warning: %d levels lack a colormap match; assigning %s. Missing levels: %s",
      length(missing_levels), missing_colour, paste(missing_levels, collapse = ", ")
    ))
  }
  add_palette <- stats::setNames(rep(missing_colour, length(missing_levels)), missing_levels)
  palette <- c(base_palette, add_palette)

  # Include the explicit "(missing)" label only if present
  has_missing_label <- n_missing_vals > 0
  if (!has_missing_label) {
    palette <- palette[names(palette) != "(missing)"]
    used_levels <- setdiff(used_levels, "(missing)")
  } else if (!"(missing)" %in% names(palette)) {
    palette <- c(palette, "(missing)" = missing_colour)
  }

  # Keep only colours for used levels and order legend labels.
  # If labels (excluding "(missing)") are numeric strings, sort numerically; otherwise alphabetical.
  used_levels <- intersect(used_levels, names(palette))
  has_missing_label <- any(used_levels == "(missing)")
  core_levels <- setdiff(used_levels, "(missing)")
  suppressWarnings({ core_nums <- as.numeric(core_levels) })
  if (length(core_levels) > 0 && all(!is.na(core_nums))) {
    ord <- order(core_nums, na.last = NA)
    ordered_core <- core_levels[ord]
  } else {
    ordered_core <- sort(core_levels)
  }
  ordered_levels <- if (has_missing_label) c(ordered_core, "(missing)") else ordered_core
  palette <- palette[ordered_levels]

  joined_dt[, (color_column) := factor(get(factor_col), levels = ordered_levels)]

  list(data = joined_dt, color_values = palette, color_column = color_column)
}

theme_00 <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(size = rel(1), margin = margin(0,0,5,0), hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.70), hjust = 0.5),
      #axis
      axis.title = element_text(size = rel(0.85)),
      axis.text = element_text(size = rel(0.70)),
      #legend
      legend.title = element_text(size = rel(0.85)),
      legend.text = element_text(size = rel(0.70)),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(1.5, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      #text
      #text=element_text(family=font, colour = "black")
    )
}


auto_plot_dimension <- function(x, y, base_dim = 4, scale_factor = 0.1, min_dim = 4, max_dim = 15) {
  range_x <- diff(range(x, na.rm = TRUE))
  range_y <- diff(range(y, na.rm = TRUE))
  max_span <- max(range_x, range_y)
  if (!is.finite(max_span) || max_span <= 0) {
    max_span <- 1
  }
  plot_dim <- base_dim + max_span * scale_factor
  min(max(plot_dim, min_dim), max_dim)
}

draw_umap_plot <- function(
  plot_data,
  color_values,
  color_column,
  subtitle = NULL,
  x_col = "UMAP1",
  y_col = "UMAP2"
) {
  p <- ggplot(plot_data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(color = !!sym(color_column)), size = 0.25, alpha = 0.35, shape = 16) +
    labs(title = "UMAP", subtitle = subtitle, x = x_col, y = y_col) +
    coord_equal() +
    theme_00(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.width = grid::unit(0.6, "lines"),
      legend.key.height = grid::unit(0.5, "lines")
    ) +
    scale_color_manual(values = color_values, na.value = "grey70", drop = TRUE) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3, alpha = 1)))

  return(p)
}
