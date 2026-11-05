# Colorectal adenoma-carcinoma progression and polyphenol prioritization

Hiu Mai

**[Read the Jupyter notebook](https://github.com/hiuuuuu33/colorectal-polyphenol-research/blob/main/Analysis.ipynb) | [Open the HTML version](https://hiuuuuu33.github.io/colorectal-polyphenol-research/)**

I wanted to understand which gene-expression changes persist as normal colorectal tissue progresses through adenoma to early cancer, and whether those changes could help prioritize polyphenols for further study. I analyzed the three tissue groups, identified genes that changed in the same direction across both transitions, and compared their positions in the human protein interaction network with recorded polyphenol targets.

This repository accompanies my manuscript, *Network-based discovery of bioactive polyphenols based on transcriptomic changes across the colorectal adenoma-carcinoma sequence*. I have included the data, my original analysis files and figures, together with an executable notebook. The notebook follows the analysis one step at a time: a short explanation, the code, then its output.

## Research summary

The expression results point to changes that begin at the adenoma stage: inflammatory recruitment increases while differentiated epithelial and metabolic functions decline. The carcinoma comparison adds extracellular-matrix organization, focal adhesion and cell-cycle activity. I identified 29 genes that change in the same direction across both sequential comparisons and used their protein relationships to connect these findings to a compound screen.

Eleven unique compounds met the network-proximity cutoff. Kaempferol ranked first, followed by 4-methylcatechol, quercetin, punicalagin and umbelliferone. I took the ten highest-ranked compounds into structural assessment. Kaempferol also had the most favorable ADH1C docking score, making it my leading candidate for experimental follow-up. The [notebook's research abstract](https://hiuuuuu33.github.io/colorectal-polyphenol-research/#research-abstract) expands this summary, and the explanations beside each figure connect the analysis to the findings in my manuscript.

| Analysis | Result |
|:--|:--|
| Expression dataset | GSE228012; 33 samples, 11 per tissue group |
| Differential expression | 725 adenoma-normal; 312 cancer-adenoma; 1,544 cancer-normal |
| Progression-associated genes | 29 genes: 8 up/up and 21 down/down |
| Compound screen | 66 recorded names representing 65 PubChem identities; 11 unique identities below the original Z-score cutoff |
| Highest-ranked candidate | Kaempferol |
| Docking records | 10 primary compounds and one supplementary structure against two receptor systems |

In the notebook, I recompute expression results from the processed FPKM matrix, protein centrality from the archived network, and observed distances between disease proteins and compound targets. For enrichment plots, I use the complete original GO and KEGG result tables. The historical compound ranking uses the saved random-distribution moments. I display the original STRING PPI export and native STRING images for each compound network. For docking, I rebuild the score table from the original workbook and show the four original images from my docking presentation individually.

## Why I chose these docking systems

I selected CXCL8 and ADH1C because both belong to the progression-associated gene set, have relatively high connectivity in its PPI network (13 and seven connections), and have experimentally determined human protein structures available. I used [6WZM](https://www.rcsb.org/structure/6WZM) and [1HT0](https://www.rcsb.org/structure/1HT0), resolved by X-ray diffraction at 2.28 and 2.00 angstrom, respectively. The notebook records their deposited refinement statistics, chain composition and receptor preparation, alongside the selection rationale in Section 3.5 of my manuscript.

All ten prioritized compounds were docked against both receptor systems. I show phloretin-6WZM and kaempferol-1HT0 as representative pairs because each had the most favorable Vina score within its receptor: -6.1 and -8.3 kcal/mol, respectively. Kaempferol also ranked first in network proximity. The [docking section](https://hiuuuuu33.github.io/colorectal-polyphenol-research/#docking) explains these choices before presenting the score table and original images.

## Files to start with

- [Analysis.ipynb](Analysis.ipynb): the Jupyter notebook, with code and saved outputs.
- [HTML notebook](https://hiuuuuu33.github.io/colorectal-polyphenol-research/): a browser-readable copy of the analysis.
- [Figure index](docs/figure-index.md): the inputs and scripts behind each manuscript figure.
- [Data and code provenance](docs/provenance.md): where the files came from and how I reconstructed the archived analyses.

```text
Analysis.ipynb          Jupyter notebook with explanations, code and outputs
data/                   Expression, enrichment, network, chemical and docking inputs
scripts/R/              Expression, enrichment, network and docking modules
scripts/python/         Network analysis and fresh randomization
scripts/                Rendering, dependencies and figure preparation
results/                Result tables, numerical checks and environment records
figures/submitted/      Nine figures extracted from my final PDF
figures/original/       Original plots and editable figure sources
figures/reproduced/     Figures generated by the analysis modules
figures/notebook/       Additional exported plots
archive/                Original R scripts, exported tables and R workspaces
upstream/               Third-party network code and its license
environment/            Dependencies and tested package versions
metadata/               File manifests and SHA-256 checksums
docs/                   HTML report and supporting documentation
report/                 Notebook stylesheet
```

I kept earlier exports in `archive/` and selected the inputs used by the final analysis under `data/`, including the 66-entry network output. Duplicate ZIP files, downloaded reference articles, system files, unrelated personal documents and intermediate manuscript drafts remain outside this public repository.

## Run the analysis

The notebook uses an R kernel, with R 4.4 or later and Python 3.12 or later. Tested package versions are listed in [environment/](environment/).

From the repository root, install the dependencies and open the notebook:

```bash
Rscript scripts/install_dependencies.R
python3 -m pip install -r environment/requirements-notebook.txt
Rscript scripts/install_notebook_kernel.R
python3 -m jupyterlab Analysis.ipynb
```

To execute every cell and export the HTML version:

```bash
python3 scripts/render_notebook.py
```

This updates the saved outputs in `Analysis.ipynb` and writes `docs/index.html`. The HTML contains its images and can be read offline. I included the input files, cached PubChem depictions and STRING responses, so the default build runs without querying those databases. Set `PYTHON` to the Python executable you want the R code to use.

The analysis modules can also run separately:

```bash
Rscript scripts/R/01_transcriptomics.R
Rscript scripts/R/02_enrichment_figures.R
python3 scripts/python/network_analysis.py --mode validate
python3 scripts/python/fetch_string_compounds.py --offline
Rscript scripts/R/03_network.R
Rscript scripts/R/04_docking.R
python3 scripts/build_chemical_structures.py
```

To generate a fresh seeded null distribution:

```bash
python3 scripts/python/network_analysis.py --mode randomize --n-random 1000
```

That command saves the new null vectors, seed and run metadata separately from the historical ranking. The docking archive includes PDB structures, search grids, scores and visualizations. A new Vina search also requires the prepared PDBQT files and output pose coordinates listed in [the docking record](docs/docking-provenance.md).

## Sources

- Expression data: [GEO GSE228012](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228012), associated with [Lin et al.](https://doi.org/10.1002/advs.202410058).
- Network framework and compound-target resource: [do Valle et al. (2021)](https://doi.org/10.1038/s43016-021-00243-7), with [their original code](https://github.com/italodovalle/polyphenols).
- Protein associations: [STRING](https://string-db.org/), using archived progression-network edges and cached version-12.0 compound-network responses.
- Chemical identifiers and depictions: [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
- Structural records: [PDB 1HT0](https://www.rcsb.org/structure/1HT0) and [PDB 6WZM](https://www.rcsb.org/structure/6WZM).

I retained the upstream polyphenol code's MIT license and attribution. Third-party data and structural records retain their source terms, as listed in [the third-party notices](THIRD_PARTY_NOTICES.md).

I organized this repository from my existing project files in November 2026, with Codex assistance for the reproduction scripts, notebook and documentation. My original files are preserved separately from those additions. My submitted manuscript describes my research contributions, mentor guidance and use of computational assistance; [the provenance record](docs/provenance.md) describes the repository preparation.
