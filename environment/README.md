# Tested Environment

The executable report was tested with R 4.4.3 and Python 3.12.14 on macOS arm64. Exact directly used R package versions are recorded in `r-packages.csv`; full R environment records are saved with each result module. Python versions are pinned in `requirements.txt`.

Pandoc is required for the HTML render. It is bundled with RStudio. A separate Pandoc installation can be selected through `RSTUDIO_PANDOC`, pointing to the folder containing the `pandoc` executable. Python can be selected through `PYTHON`.

The default build runs exclusively on included input files and cached source responses. Package installation is a separate step that requires internet access. Figures 1-6 and the archived docking-score summary do not require any remote database queries. The full 1,000-iteration network randomization is an optional new computation and takes longer than the default historical reconstruction.

For PDF figure extraction only, install PyMuPDF separately and run `python3 scripts/extract_submitted_figures.py`. The nine extracted figures are already included in the repository.
