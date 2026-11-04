# Install the R packages required to run the analysis and render the report.

cran <- c("rmarkdown", "knitr", "ggplot2", "ggrepel", "pheatmap", "patchwork", "readxl", "jsonlite", "digest")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
message("R dependencies ready. Pandoc must also be available (bundled with RStudio or installed separately).")
