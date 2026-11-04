# Tested Environment

The executable report was tested with R 4.4.3 and Python 3.12.14 on macOS arm64. Exact directly used R package versions are recorded in `r-packages.csv`; full R environment records are saved with each result module. Python versions are pinned in `requirements.txt`.

The Jupyter HTML export uses nbconvert and does not require Pandoc. The alternative R Markdown render requires Pandoc, which is bundled with RStudio. A separate Pandoc installation can be selected through `RSTUDIO_PANDOC`, pointing to the folder containing the `pandoc` executable. Python can be selected through `PYTHON`.

The default build runs exclusively on included input files and cached source responses. Package installation is a separate step that requires internet access. Figures 1-6 and the archived docking-score summary do not require any remote database queries. The full 1,000-iteration network randomization is an optional new computation and takes longer than the default historical reconstruction.

For PDF figure extraction only, install PyMuPDF separately and run `python3 scripts/extract_submitted_figures.py`. The nine extracted figures are already included in the repository.

## Jupyter Notebook

`Analysis.ipynb` runs with the R kernel, so the expression model and plotting code use the same R packages as the R Markdown report. The saved notebook includes its executed outputs and can be read on GitHub without installing Jupyter.

To run it locally from the repository root:

```bash
python3 -m pip install -r environment/requirements-notebook.txt
Rscript scripts/install_dependencies.R
Rscript scripts/install_notebook_kernel.R
python3 -m jupyterlab Analysis.ipynb
```

The kernel is registered for the current user as `crc-polyphenol-r`, displayed as `R (colorectal polyphenols)`. Set `PYTHON` to the intended Python executable before installing the kernel and running the notebook if several Python environments are available. This leaves other registered R kernels available.

The notebook environment was tested with IRkernel 1.3.2 and IRdisplay 1.1. The Python notebook packages are pinned separately in `requirements-notebook.txt`; the research packages remain in `requirements.txt`.

Run `python3 scripts/render_notebook.py` to execute all cells and update the saved notebook and `docs/index.html`. The exporter embeds the output images and styles and disables external script dependencies, so the HTML can also be opened offline.
