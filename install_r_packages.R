# Install CartLoader R script dependencies
#
# Run with: `Rscript cartloader/r/install_r_dependencies.R`

required_packages <- c(
  "argparse",
  "data.table",
  "RcppParallel",
  "ggplot2",
  "uwot"
)

installed <- rownames(installed.packages())
missing <- setdiff(required_packages, installed)

if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing)
} else {
  message("All required packages are already installed.")
}
