args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
root <- normalizePath(root, mustWork = TRUE)
if (!file.exists(file.path(root, "report", "analysis.Rmd"))) stop("Run from the repository root.")
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Run Rscript scripts/install_dependencies.R first.")
if (!rmarkdown::pandoc_available()) {
  stop("Pandoc is required. Install Pandoc or render using RStudio, which includes it.")
}
dir.create(file.path(root, "docs"), showWarnings = FALSE)
Sys.setenv(RESEARCH_REPO_ROOT = root)
rmarkdown::render(file.path(root, "report", "analysis.Rmd"),
                  output_file = "index.html", output_dir = file.path(root, "docs"),
                  knit_root_dir = root, envir = new.env(parent = globalenv()),
                  clean = TRUE, quiet = FALSE)
file.create(file.path(root, "docs", ".nojekyll"))
message("Report ready: ", file.path(root, "docs", "index.html"))
