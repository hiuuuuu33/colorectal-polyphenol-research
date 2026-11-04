# Network inputs and analysis lineage

The files in this directory are copies. Original project files were left in place.

| File | Original project location | Role |
| --- | --- | --- |
| `HumanInteractome_v2017.csv` | `polyphenols-master/data/` | Human protein interactome distributed with the do Valle software |
| `PolyphenolProteinInteractions.csv` | `polyphenols-master/data/` | Compound-target records, with PubChem and Entrez identifiers |
| `CRA_diseasegenes.csv` | `polyphenols-master/data/` | Disease-gene input used for the historical proximity run |
| `carcinogenesis_with_entrez.csv` | `polyphenols-master/data/` | Original symbol-to-Entrez mapping intermediate |
| `proximity_original_66_records.csv` | `polyphenols-master/out/merged_zscore_proximity.csv` | Historical 66-record result used by the submitted manuscript |
| `ppi_original_cytoscape_nodes.csv` | `数据分析/nodecsv.csv` | Historical Cytoscape centrality export |
| `string_node_degrees.tsv` | `数据分析/网络/` | All 29 STRING query nodes, including isolated nodes |
| `string_interactions_short.tsv` | `数据分析/网络/` | Saved 61-edge progression-associated STRING network |

The 42-record output stored separately in the old `数据分析` directory was not used
for the submitted result. Resveratrol and trans-resveratrol share PubChem CID
445154 and are merged for the unique-identity ranking. The Figure 6A distribution
retains all 66 historical records; the Figure 6B ranking uses 65 identities.

The archived disease input contains 28 identifiers, of which 24 map into the
interactome. The validation script intentionally reconstructs this historical
input so its distance values can be checked against the submitted analysis.
It preserves the interactome's self-loops because they affect degree-bin
matching. It does not silently change these inputs.

## Upstream software and data credit

Original software: [italodovalle/polyphenols](https://github.com/italodovalle/polyphenols).
The bundled upstream code retains its MIT license, copyright (c) 2020 italodovalle,
in `upstream/polyphenols/LICENSE`. The local upstream snapshot has no Git commit
metadata; its exact remote commit is unrecorded. File checksums document the
actual local snapshot included here.

Citation: do Valle, I. F. et al. Network medicine framework shows that proximity
of polyphenol targets and disease proteins predicts therapeutic effects of
polyphenols. *Nature Food* 2, 143-155 (2021).

The upstream README identifies STITCH experimental interactions as the source
of the compound-target data. The individual source database terms and citations
continue to apply to those data. The MIT notice is retained for the software;
no new license claim is made over third-party database contents.

## Reproduction modes

`python scripts/python/network_analysis.py` independently recomputes all 66
observed closest and mean shortest-path distances, reconstructs PPI centralities,
and checks them against the original exports. It also recalculates each Z score
as `(recalculated_closest - saved_null_mean) / saved_null_standard_deviation`
and checks it against the saved Z. `Rscript scripts/R/03_network.R`
redraws figures from the saved and recalculated tables.

The historical run saved null means and standard deviations but did not retain
individual randomizations. Therefore, its Z values are plotted from the archived
CSV. A fresh optional run is available:

```sh
python scripts/python/network_analysis.py --mode randomize --n-random 1000 --seed 452456
```

This writes every new closest-distance null value and run metadata under
`results/network/fresh_randomization/`. It uses the upstream degree-bin scheme
with explicitly sorted node inputs, preserving self-loops and the separate
same-seed sampling streams. New randomizations are stored separately from the
historical scores. They can differ from the old result because the historical
software environment and node ordering were not fully recorded.

## Compound networks (Figure 8)

`scripts/python/fetch_string_compounds.py` retrieves the ten candidate target
networks from the version-fixed STRING 12.0 API. The target list comes from the
original compound-target CSV, not just targets mapped to the proximity network.
Parameters follow the submitted Figure 8 caption: human, functional associations,
combined score at least 0.400, no added proteins. A singleton remains a single
node. Frozen identifier mappings, interaction tables and request metadata support
offline replotting. These are newly reconstructed networks; the submitted PDF
preserves the original figure appearance.

API documentation: <https://version-12-0.string-db.org/help/api/>.
