# Install and register the R kernel used by Analysis.ipynb.

python <- Sys.getenv("PYTHON", unset = Sys.which("python3"))
if (!nzchar(python)) stop("Python 3 is required; set PYTHON to its executable.")
python <- normalizePath(python, mustWork = TRUE)
Sys.setenv(PATH = paste(dirname(python), Sys.getenv("PATH"), sep = .Platform$path.sep))
status <- system2(python, c("-m", "jupyter", "kernelspec", "--version"))
if (status != 0) stop("Install environment/requirements-notebook.txt in the selected Python environment first.")
if (!requireNamespace("IRkernel", quietly = TRUE)) {
  install.packages("IRkernel", repos = "https://cloud.r-project.org")
}
status <- IRkernel::installspec(user = TRUE, name = "crc-polyphenol-r",
                              displayname = "R (colorectal polyphenols)")
if (status != 0) stop("The R kernel could not be registered.")
message("Notebook kernel ready: crc-polyphenol-r")
