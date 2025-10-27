# Install CartLoader R script dependencies
#
# Run with:
#   Rscript install_r_packages.R
#   Rscript install_r_packages.R --libdir /home/weiqiuc/R/x86_64-pc-linux-gnu-library/4.4
#
# Supported args:
#   --libdir PATH   Install/lookup packages in PATH (created if missing)
#
# Policy: install only strong deps (Depends/Imports/LinkingTo). Skip Suggests & vignettes.

# Use a CRAN mirror explicitly to avoid interactive prompts
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Avoid heavy extras
options(build_vignettes = FALSE, Ncpus = max(1L, parallel::detectCores(logical = TRUE)))
Sys.setenv(R_BUILD_VIGNETTES = "0")

# --- simple arg parsing (no external deps) ---
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  if (length(args) == 0) return(default)
  eq_idx <- grep(paste0("^", flag, "="), args)  # --flag=value
  if (length(eq_idx) == 1L) return(sub(paste0("^", flag, "="), "", args[eq_idx]))
  pos_idx <- which(args == flag)                # --flag value
  if (length(pos_idx) == 1L && pos_idx < length(args)) return(args[pos_idx + 1L])
  default
}

# Determine target user library
user_lib <- get_arg("--libdir", default = NA_character_)
if (is.na(user_lib) || user_lib == "") {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (user_lib == "") {
    user_lib <- file.path(
      Sys.getenv("HOME"), "R",
      paste0(R.version$platform, "-", getRversion()[, "major"], ".", getRversion()[, "minor"])
    )
  }
} else {
  user_lib <- path.expand(user_lib)
}
user_lib <- path.expand(user_lib)
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)

# Prepend user library so installs go there
.libPaths(c(user_lib, .libPaths()))
message("Using library: ", .libPaths()[1])

required_packages <- c("argparse", "data.table", "RcppParallel", "ggplot2", "uwot")

# Consider all libraries in .libPaths() when checking what is installed
installed <- rownames(installed.packages(lib.loc = .libPaths()))
missing <- setdiff(required_packages, installed)

if (length(missing) > 0) {
  message("Installing missing packages (strong deps only): ", paste(missing, collapse = ", "))
  # Guard: ensure the target lib is writable
  if (file.access(.libPaths()[1], 2) != 0) {
    stop("First library path '", .libPaths()[1], "' is not writable. ",
         "Use --libdir to choose a writable directory or set R_LIBS_USER.")
  }
  # Only strong deps (Depends/Imports/LinkingTo); skip Suggests
  install.packages(
    missing,
    lib = .libPaths()[1],
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
} else {
  message("All required packages are already installed.")
}