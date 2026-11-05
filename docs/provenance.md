# Data, Code and Reconstruction Provenance

## Source Preservation

All material was copied into a separate repository. The original research directory was read without moving, renaming, replacing or deleting its files. A complete private source inventory and a before/after checksum comparison are retained outside this repository. The public file manifests use repository-relative paths.

The manuscript is kept outside this repository. `metadata/submitted_figures.csv` records the embedded figures extracted from it. `metadata/original_archive_manifest.csv` connects preserved historical scripts, tables, workspaces and figure sources to their original relative paths.

## Original Research Files

- Original R scripts are in `archive/original-r/`, including historical variants.
- Original exported tables are in `archive/original-tables/`; some belong to earlier analysis versions. They remain labeled as archived exports.
- Original R workspaces are in `archive/r-snapshots/`.
- Original generated images and editable figure sources are in `figures/original/`.
- Original molecular-docking records are in `data/docking/`.
- Original polyphenol code and license are retained in `upstream/`.

The current workflow reads explicitly selected files from `data/`. In particular, the final network analysis uses `polyphenols-master/out/merged_zscore_proximity.csv`, containing 66 entries. The earlier 42-entry export is retained as a historical table in the archive and is outside the executable analysis.

## New Reproduction Work

The modules under `scripts/R/01_*` to `04_*`, the Python network reconstruction, the chemical depiction assembly, the Jupyter notebook and repository documentation were assembled in November 2026 with Codex assistance. Their purpose is to make the supplied research executable and traceable. The analysis record is a retrospective reconstruction of the archived workflow. File modification timestamps have not been converted into a claimed day-by-day research diary.

`Analysis.ipynb` uses a native R kernel. It writes out the expression preprocessing, limma model and progression-gene intersection, followed by plotting cells and calls to the archived-result modules. Its execution counts, tables and PNG outputs come from a complete kernel execution. `scripts/render_notebook.py` reruns the cells before exporting the same inputs and outputs to `docs/index.html`.

1. Expression analysis starts from the processed FPKM matrix and reproduces the historical unblocked limma model. All three differential-expression exports and the 29-gene identities are checked numerically.
2. Functional-enrichment figures use 15 complete frozen result tables. This preserves their historical annotation state; the default render filters and plots those results.
3. Protein centrality and observed proximity distances are recomputed from the supplied graphs and mappings. The historical mapping contains 24 disease nodes; the expression-derived set contains 29 genes. Historical null-distribution moments are retained with the source output. A separate fresh-randomization mode saves newly generated null vectors and metadata.
4. Chemical depictions are fetched by PubChem CID once and then cached. The progression PPI displays the original STRING export, checked against all 29 nodes and 61 archived edges. Compound target networks use STRING 12.0, human taxonomy 9606 and combined score 0.400, with no added neighbors. Native STRING images are cached with the responses and request parameters, and their node and edge sets are checked against the analysis tables. Each compound network is displayed separately.
5. Docking score analyses read the original workbook. The four original PNG panels in `data/docking/records/docking_results.pptx` are extracted without modifying their bytes and displayed individually in the notebook; their sources and hashes are recorded in `metadata/docking_panels.csv`. Source PDB structures and grids are retained. Prepared receptor/ligand PDBQT and output pose coordinates were absent from the supplied archive; repeating the docking search requires those additional inputs.

The final PDF is the source of the submission's wording. This repository documents the input files that reproduce the archived outputs without modifying the submitted manuscript.

## Excluded Material

Downloaded full-text reference articles are cited by their source URLs rather than republished here. Duplicate ZIP archives, unrelated personal documents, Office lock files, operating-system metadata and intermediate manuscript drafts remain in the original folder. Connectivity-map example datasets unused by the colorectal workflow remain with the original upstream download.

## Remaining Inputs

The complete historical randomization vectors and prepared docking input/output coordinates were absent from the supplied source directory. Newly generated null vectors are saved independently from the historical ranking. Frozen enrichment results preserve the available record of the original annotation calculations.
