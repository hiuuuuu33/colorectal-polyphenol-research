# Molecular docking: records and reproducible outputs

The submitted manuscript reports the docking study in Section 3.5, Table 5 and
Figure 9. This repository preserves the original structural inputs and result
workbook and regenerates the score table and an additional score-summary plot.
The submitted manuscript itself is unchanged.

## Figure and table map

| Submitted item | Contents | Repository treatment |
| --- | --- | --- |
| Table 5 | Best-mode Vina scores for ten prioritized compounds and one supplementary punicalagin structure | Rebuilt directly from `data/docking/records/docking_results.xlsx` |
| Figure 9A-B | Phloretin in the 6WZM CXCL8-Fab assembly; overall view and two-dimensional interaction map | Original PNG panels from slide 1, displayed separately in the notebook |
| Figure 9C-D | Kaempferol in the 1HT0 ADH1C model; overall view and two-dimensional interaction map | Original PNG panels from slide 2, displayed separately in the notebook |
| Additional repository plot | All 22 archived scores, separated by receptor system | Generated as `figures/reproduced/docking/docking_score_summary.png` and PDF |

Figure 8 in the final submission contains the compound-associated target
networks, rather than docking results.

## Included source records

All source files below were copied from the author's original project directory;
the originals remain untouched.

| Repository path | Source item |
| --- | --- |
| `data/docking/raw_structures/ADH1C_1HT0.pdb` | Original RCSB structure download, PDB [1HT0](https://www.rcsb.org/structure/1HT0) |
| `data/docking/raw_structures/CXCL8_6WZM.pdb` | Original RCSB structure download, PDB [6WZM](https://www.rcsb.org/structure/6WZM) |
| `data/docking/grids/ADH1C_1HT0.txt` | Archived 1HT0 search-box parameters |
| `data/docking/grids/CXCL8_6WZM.txt` | Archived 6WZM search-box parameters |
| `data/docking/records/docking_results.xlsx` | Original score workbook, Sheet1; Sheets2-3 are empty |
| `data/docking/records/docking_results.pptx` | Archived visualization slides |
| `data/docking/original_visualizations/*.png` | Four archived visualization exports, named by receptor, PDB identifier and view |
| `figures/original/docking/*.png` | Four original PNG media files extracted byte for byte from the archived PPT |
| `metadata/docking_panels.csv` | Slide and media mapping, original dimensions and SHA-256 digests for the displayed panels |

The notebook displays the four full-resolution PPT images corresponding to
submitted Figure 9. Each image retains its original pixels and white background;
the images have not been redrawn, cropped or rescaled on disk. The separate
exports in `data/docking/original_visualizations/` include earlier views. Their
filenames identify the receptor, PDB entry and 3D or 2D view. The original names
are recorded in `metadata/original_archive_manifest.csv`; all image bytes are
preserved.

[Phloretin / 6WZM: 3D image](../figures/original/docking/phloretin_6wzm_3d.png) |
[Interaction map](../figures/original/docking/phloretin_6wzm_interactions.png)

[Kaempferol / 1HT0: 3D image](../figures/original/docking/kaempferol_1ht0_3d.png) |
[Interaction map](../figures/original/docking/kaempferol_1ht0_interactions.png)

To extract the original panels again:

```sh
python3 scripts/extract_docking_panels.py
```

## Why these targets and structures

I selected CXCL8 and ADH1C for structural assessment because they occur in the
29-gene progression-associated set, are relatively well connected in its PPI
network, and have experimentally determined human protein structures available.
CXCL8 has 13 direct connections in this network and ADH1C has seven. These are
the selection criteria reported in Section 3.5 of my submitted manuscript.

The deposited structures provide experimental atomic coordinates for receptor
preparation. Their quality is documented by the resolution and refinement
statistics below, checked against the preserved PDB headers and the linked RCSB
entries. R factors are the depositor-reported values.

| Structural feature | CXCL8-Fab / [6WZM](https://www.rcsb.org/structure/6WZM) | ADH1C / [1HT0](https://www.rcsb.org/structure/1HT0) |
| --- | --- | --- |
| Experimental method | X-ray diffraction | X-ray diffraction |
| Resolution | 2.28 angstrom | 2.00 angstrom |
| R-work | 0.187 | 0.193 |
| R-free | 0.229 | 0.231 |
| Human protein | CXCL8 | ADH1C gamma-2 |
| Deposited protein chains | CXCL8 E/F; Fab heavy A/C, light B/D | ADH1C A/B |

All ten prioritized compounds were docked against both receptor systems. The
representative Figure 9 pairs were chosen from those results: phloretin has the
most favorable primary-compound score in 6WZM (-6.1 kcal/mol), and kaempferol
has the most favorable score in 1HT0 (-8.3 kcal/mol). Kaempferol also leads the
network-proximity ranking. This explains the choice of the two illustrated
poses separately from the choice of the two receptors.

## Structural systems and recorded settings

The original 1HT0 file is human ADH1C gamma-2 and contains two protein chains,
two NAD molecules and four zinc ions. The submitted Methods report removal of
NAD and zinc during preparation. The search-box centre coincides with the
crystallographic NAD A1377 site.

The original 6WZM file contains a CXCL8-Fab assembly. Chains A/C are Fab heavy
chains, B/D are Fab light chains and E/F are CXCL8. The submitted Methods report
using the full chain A-F assembly. Its search-box centre coincides with the
crystallographic tetraethylene glycol site, PG4 F101, at the chain C/F interface.

| Parameter | ADH1C / 1HT0 | CXCL8-Fab / 6WZM |
| --- | ---: | ---: |
| `center_x` (angstrom) | 20.446 | 23.673 |
| `center_y` (angstrom) | -4.700 | -1.174 |
| `center_z` (angstrom) | 4.799 | 32.063 |
| `size_x` (angstrom) | 27.580 | 14.306 |
| `size_y` (angstrom) | 27.580 | 14.306 |
| `size_z` (angstrom) | 27.580 | 14.306 |
| `exhaustiveness` | 8 | 8 |
| `num_modes` | 9 | 9 |
| `energy_range` (kcal/mol) | 3 | 3 |

The submitted Methods specify AutoDock Tools 1.5.6 for receptor/ligand
preparation, Open Babel for ligand conversion, AutoDock Vina 1.1.2 for docking,
and BIOVIA Discovery Studio Visualizer for interaction visualization. The
Visualizer version and random seed are not recorded in the available files.
The raw PDBs here preserve the downloaded structures, before the reported
preparation steps.

## Run the reproducible record analysis

From the repository root:

```sh
Rscript scripts/R/04_docking.R
```

Or within R or the Jupyter notebook:

```r
source("scripts/R/04_docking.R")
docking <- run_docking(".")
docking$scores
docking$plot
```

R dependencies: `readxl`, `ggplot2`. The module reads all 11 workbook rows and
extracts PubChem CIDs from their archived pose labels. It verifies record
identities, orders the ten primary records as in Table 5, places the additional
punicalagin CID 44584733 last, and writes:

- `results/docking/table5_docking_scores.csv`
- `results/docking/archived_grid_parameters.csv`
- `results/docking/reproduction_status.csv`
- `results/docking/sessionInfo.txt`
- `figures/reproduced/docking/docking_score_summary.png`
- `figures/reproduced/docking/docking_score_summary.pdf`

These outputs reproduce the archived score analysis. Generating new Vina scores
and interaction diagrams requires the prepared receptor and ligand PDBQT files,
the final pose coordinates, preparation settings, commands and seed. Those files
were not present in the provided project directory. The workbook's `.pdbqt`
strings are record labels; the underlying coordinate files are not included.
Once recovered, they can be added alongside the existing raw structures without
changing the archived score records.

## Score interpretation

The table contains ten primary compound records plus a supplementary
stereospecified punicalagin record, rather than eleven independent prioritized
compounds. Scores are preserved exactly as recorded, including the positive
6WZM score for the supplementary structure. Comparisons are made within each
prepared receptor system. The additional summary plot is a new reproducible
display of the existing scores, not a replacement for the submitted pose figure.

## References

- RCSB Protein Data Bank: [1HT0](https://www.rcsb.org/structure/1HT0),
  [6WZM](https://www.rcsb.org/structure/6WZM).
- Trott O, Olson AJ. AutoDock Vina. *Journal of Computational Chemistry*.
  2010;31:455-461. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334).
