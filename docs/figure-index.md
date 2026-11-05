# Manuscript Figure and Table Index

The reference is the unchanged 29-page `manuscript/submitted-manuscript.pdf`. All nine embedded figures are preserved in `figures/submitted/`.

The [HTML notebook](https://hiuuuuu33.github.io/colorectal-polyphenol-research/) displays the executable code immediately before each figure. Its source is [Analysis.ipynb](../Analysis.ipynb) at the repository root, with executed R cells and embedded PNG outputs.

| Submitted item | PDF page | Inputs | Executable companion | Reproduction path |
|:--|--:|:--|:--|:--|
| Figure 1: PCA and volcano plots | 8 | `data/transcriptomics/` | `scripts/R/01_transcriptomics.R` | Recomputed from FPKM |
| Figure 2: progression-associated genes | 9 | Sequential expression contrasts | `scripts/R/01_transcriptomics.R` | Recomputed 29-gene intersection, heatmap and profiles |
| Figure 3: adenoma-normal enrichment | 10 | `data/enrichment/frozen/` | `scripts/R/02_enrichment_figures.R` | Regenerated from complete frozen enrichment results |
| Figure 4: cancer-adenoma enrichment | 12 | `data/enrichment/frozen/` | `scripts/R/02_enrichment_figures.R` | Regenerated from complete frozen enrichment results |
| Figure 5: progression PPI topology | 13 | Archived STRING edge and node tables; original STRING SVG and PNG | `scripts/python/network_analysis.py`, `scripts/build_string_visuals.py` | Original native STRING image on white; its 29 nodes and 61 edges checked against the analysis tables; recomputed centrality |
| Figure 6: compound proximity | 15 | Archived interactome, targets and null moments | Same network scripts | Observed distances and Z-score arithmetic checked; regenerated ranking |
| Figure 7: chemical structures | 18 | Top-ten PubChem CIDs; cached PNG depictions | `scripts/build_chemical_structures.py` | Rebuilt companion from cached source depictions |
| Figure 8: compound target networks | 19 | Original target sets; cached STRING 12.0 mapping, edges and native images | `scripts/build_string_visuals.py` and frozen query responses | Ten separate native STRING images; image node and edge sets checked against the cached tables |
| Table 5: docking scores | 20 | Original docking XLSX | `scripts/R/04_docking.R` | Rebuilt score table and summary plot |
| Figure 9: representative poses | 21 | Four original PNG panels embedded in the docking PPT | `scripts/extract_docking_panels.py` | Original panel bytes extracted without cropping or redrawing; full-size 3D and 2D images displayed separately |

The regenerated companions use readable report layouts. Original figure compositions remain available alongside them. New compound-network responses retain their access dates and parameters in `data/network/string_visuals/`. `results/network/string_visuals.json` records the image paths and graph checks. Compact notebook previews link to the full-resolution native PNGs. `metadata/docking_panels.csv` records the PPT media sources, image sizes and SHA-256 digests of the four original docking panels.
