# Manuscript Figure and Table Index

The reference is the unchanged 29-page `manuscript/submitted-manuscript.pdf`. All nine embedded figures are preserved in `figures/submitted/`.

| Submitted item | PDF page | Inputs | Executable companion | Reproduction path |
|:--|--:|:--|:--|:--|
| Figure 1: PCA and volcano plots | 8 | `data/transcriptomics/` | `scripts/R/01_transcriptomics.R` | Recomputed from FPKM |
| Figure 2: progression-associated genes | 9 | Sequential expression contrasts | `scripts/R/01_transcriptomics.R` | Recomputed 29-gene intersection, heatmap and profiles |
| Figure 3: adenoma-normal enrichment | 10 | `data/enrichment/frozen/` | `scripts/R/02_enrichment_figures.R` | Regenerated from complete frozen enrichment results |
| Figure 4: cancer-adenoma enrichment | 12 | `data/enrichment/frozen/` | `scripts/R/02_enrichment_figures.R` | Regenerated from complete frozen enrichment results |
| Figure 5: progression PPI topology | 13 | Archived STRING edge and node tables | `scripts/python/network_analysis.py`, `scripts/R/03_network.R` | Recomputed centrality; regenerated layout |
| Figure 6: compound proximity | 15 | Archived interactome, targets and null moments | Same network scripts | Observed distances and Z-score arithmetic checked; regenerated ranking |
| Figure 7: chemical structures | 18 | Top-ten PubChem CIDs; cached PNG depictions | `scripts/build_chemical_structures.py` | Rebuilt companion from cached source depictions |
| Figure 8: compound target networks | 19 | Original target sets; cached STRING 12.0 mapping and edges | Network scripts and frozen query responses | Reconstructed with submitted query parameters; retrieval metadata retained |
| Table 5: docking scores | 20 | Original docking XLSX | `scripts/R/04_docking.R` | Rebuilt score table and summary plot |
| Figure 9: representative poses | 21 | Embedded submitted figure; original docking images | Preserved illustration | Archived visualization; prepared pose coordinates needed for a fresh docking run |

The regenerated companions use readable report layouts. Original figure compositions remain available alongside them. New compound-network responses are identified as a reconstruction and retain their access dates.
