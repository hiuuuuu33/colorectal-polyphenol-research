# Network module validation

The default validation command completed with Python 3.12.14, NetworkX 3.7,
NumPy 2.3.5 and pandas 2.2.3. Replotting completed with R 4.4.3 and the installed
ggplot2 and ggrepel packages.

| Check | Result |
| --- | --- |
| Historical records / unique PubChem identities | 66 / 65 |
| Records / identities below Z = -1.96 | 12 / 11 |
| Independently recalculated closest distances | All 66 agree; maximum absolute error 2.22e-16 |
| Independently recalculated mean shortest distances | All 66 agree; maximum absolute error 4.44e-16 |
| Z from recalculated observed distance and saved null moments | All 66 agree; maximum absolute error 2.66e-15 |
| STRING progression network | 29 nodes, 61 edges |
| PPI degree, component-normalized betweenness and closeness | Agree with the original Cytoscape export |
| Interactome largest connected component | 17,651 nodes, 351,393 edges, retaining historical self-loops |
| Fresh randomization smoke test | Kaempferol, 20 iterations, seed 452456; completed |
| Full fresh 65-identity randomization | Optional command provided; not run during packaging |

The 20-iteration output tests the execution and null-vector export. Its result
is stored separately and is not used for the submitted ranking or figures.

Figure 5A is a new layout of the saved 29-node edge list. Figure 5B displays
recalculated centralities. Figure 6A plots all 66 historical records and Figure
6B merges PubChem duplicates before selecting ten identities. All new layouts
are deterministic and preserve the source tables.

Figure 8 is newly reconstructed from version-fixed STRING 12.0 requests, using
all nonduplicated recorded compound targets, human taxon 9606, a functional
association score cutoff of 0.400, and no additional proteins. The archived
submitted PDF retains the original figure appearance. Request metadata records
retrieval times and SHA-256 checksums; edge and mapping tables allow offline
replotting. Dense network labels are abbreviated visually by overlap suppression;
the complete node names remain in the machine-readable tables.

## Commands

```sh
python scripts/python/network_analysis.py
Rscript scripts/R/03_network.R

# Optional repeat of the smoke test
python scripts/python/network_analysis.py --mode randomize --n-random 20 --chemical kaempferol

# Optional full fresh randomization
python scripts/python/network_analysis.py --mode randomize --n-random 1000 --seed 452456

# Optional online refresh of frozen STRING inputs
python scripts/python/fetch_string_compounds.py --refresh
```
