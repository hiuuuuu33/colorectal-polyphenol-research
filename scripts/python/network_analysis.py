#!/usr/bin/env python3
"""Calculate protein-network centrality and compound proximity.

Validate the saved results or generate a seeded null distribution. Degree-bin
sampling follows do Valle et al. (2021), retaining self-loop degrees. Historical
Z scores use saved randomization moments.
"""

import argparse
import csv
import hashlib
import json
import math
import platform
import random
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "network"
OUT = ROOT / "results" / "network"


def read_inputs():
    edges = pd.read_csv(DATA / "HumanInteractome_v2017.csv")
    graph = nx.from_pandas_edgelist(edges, "EntrezA", "EntrezB")
    graph = graph.subgraph(max(nx.connected_components(graph), key=len)).copy()
    disease = pd.read_csv(DATA / "CRA_diseasegenes.csv")
    targets = pd.read_csv(DATA / "PolyphenolProteinInteractions.csv")
    saved = pd.read_csv(DATA / "proximity_original_66_records.csv", index_col=0)
    target_sets = {
        name: set(frame.entrez_id.dropna().astype(int)) & set(graph)
        for name, frame in targets.groupby("chemical", sort=False)
    }
    disease_set = set(disease.entrez_id.astype(int)) & set(graph)
    return graph, disease, targets, saved, target_sets, disease_set


def original_records(targets, saved):
    cids = targets.groupby("chemical").pubchem_compound_id.first().astype(int)
    records = saved.copy()
    records["pubchem_cid"] = records.chemical.map(cids)
    assert records.pubchem_cid.notna().all()
    records.pubchem_cid = records.pubchem_cid.astype(int)
    records = records.sort_values(["z_closest", "chemical"], kind="stable")
    records.to_csv(OUT / "proximity_all_66_records.csv", index=False)
    unique = records.drop_duplicates("pubchem_cid").copy()
    unique["display_name"] = unique.chemical
    unique.loc[unique.pubchem_cid == 445154, "display_name"] = "trans-resveratrol"
    unique.loc[unique.pubchem_cid == 8530, "display_name"] = "1,4-naphthoquinone"
    unique["rank"] = range(1, len(unique) + 1)
    unique.to_csv(OUT / "proximity_65_unique_identities.csv", index=False)
    unique.head(10).to_csv(OUT / "top10_candidates.csv", index=False)
    entries = []
    for row in unique.head(10).itertuples():
        subset = targets[targets.chemical == row.chemical]
        subset = subset.drop_duplicates("entrez_id")
        for target in subset.itertuples():
            entries.append({"rank": row.rank, "chemical": row.chemical,
                            "display_name": row.display_name,
                            "pubchem_cid": row.pubchem_cid,
                            "entrez_id": int(target.entrez_id),
                            "symbol": target.symbol})
    pd.DataFrame(entries).to_csv(OUT / "top10_recorded_targets.csv", index=False)
    return unique


def recalculate_observed(graph, saved, target_sets, disease_set):
    distances = {node: nx.single_source_shortest_path_length(graph, node)
                 for node in sorted(disease_set)}
    nearest = nx.multi_source_dijkstra_path_length(graph, disease_set, weight=None)
    checks = []
    for row in saved.itertuples():
        targets = sorted(target_sets[row.chemical])
        closest = np.mean([nearest[node] for node in targets])
        shortest = np.mean([distances[d][t] for d in disease_set for t in targets])
        checks.append({"chemical": row.chemical,
                       "mapped_targets_recalculated": len(targets),
                       "mapped_disease_recalculated": len(disease_set),
                       "closest_recalculated": closest,
                       "closest_saved": row.closest,
                       "closest_absolute_error": abs(closest - row.closest),
                       "z_closest_recalculated_from_saved_moments": (closest - row.avg_closest) / row.std_closest,
                       "z_closest_saved": row.z_closest,
                       "z_closest_absolute_error": abs((closest - row.avg_closest) / row.std_closest - row.z_closest),
                       "shortest_recalculated": shortest,
                       "shortest_saved": row.shortest,
                       "shortest_absolute_error": abs(shortest - row.shortest)})
    result = pd.DataFrame(checks)
    result.to_csv(OUT / "observed_distance_validation.csv", index=False)
    assert result.closest_absolute_error.max() < 1e-12
    assert result.shortest_absolute_error.max() < 1e-12
    assert result.z_closest_absolute_error.max() < 1e-12
    assert (result.mapped_targets_recalculated.values == saved.n_mapped_chemical.values).all()
    return result


def ppi_topology():
    edges = pd.read_csv(DATA / "string_interactions_short.tsv", sep="\t")
    nodes = pd.read_csv(DATA / "string_node_degrees.tsv", sep="\t")
    graph = nx.from_pandas_edgelist(edges, "#node1", "node2")
    graph.add_nodes_from(nodes["#node"])
    rows, coordinates = [], []
    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    for component_index, component in enumerate(components):
        sub = graph.subgraph(component).copy()
        between = nx.betweenness_centrality(sub, normalized=True)
        close = nx.closeness_centrality(sub, wf_improved=False)
        pos = nx.spring_layout(sub, seed=2026, iterations=300)
        for name in sorted(component):
            x, y = pos[name]
            if component_index:
                x, y = -1.2 + component_index * 0.52 + x * 0.12, -1.45 + y * 0.12
            coordinates.append({"name": name, "x": x, "y": y})
            rows.append({"name": name, "degree": graph.degree(name),
                         "betweenness": between[name], "closeness": close[name],
                         "component_size": len(component),
                         "component": component_index + 1})
    result = pd.DataFrame(rows).sort_values(["degree", "betweenness", "name"],
                                           ascending=[False, False, True])
    result.to_csv(OUT / "ppi_centrality_recalculated.csv", index=False)
    pd.DataFrame(coordinates).to_csv(OUT / "ppi_layout.csv", index=False)
    original = pd.read_csv(DATA / "ppi_original_cytoscape_nodes.csv")
    merged = result.merge(original, on="name")
    assert np.allclose(merged.degree, merged.Degree)
    assert np.allclose(merged.betweenness, merged.BetweennessCentrality)
    assert np.allclose(merged.closeness, merged.ClosenessCentrality)
    return graph, result


def degree_bins(graph, minimum):
    by_degree = defaultdict(list)
    for node, degree in graph.degree():
        by_degree[degree].append(node)
    degrees = sorted(by_degree)
    bins, index = [], 0
    while index < len(degrees):
        low = degrees[index]
        values = sorted(by_degree[low])
        while len(values) < minimum:
            index += 1
            if index == len(degrees):
                break
            values.extend(sorted(by_degree[degrees[index]]))
        if index == len(degrees):
            index -= 1
        high = degrees[index]
        index += 1
        if len(values) < minimum:
            prev_low, _, prev_values = bins[-1]
            bins[-1] = (prev_low, high, prev_values + values)
        else:
            bins.append((low, high, values))
    return bins


def random_sets(graph, selected, bins, count, seed):
    rng = random.Random(seed)
    choices = {}
    for node in sorted(selected):
        degree = graph.degree(node)
        choices[node] = next([n for n in values if n != node]
                             for low, high, values in bins if low <= degree <= high)
    samples = []
    for _ in range(count):
        sample = set()
        for node in sorted(choices):
            chosen = rng.choice(choices[node])
            for _ in range(20):
                if chosen in sample:
                    chosen = rng.choice(choices[node])
            sample.add(chosen)
        samples.append(sample)
    return samples


def fresh_randomization(graph, unique, target_sets, disease_set, count, seed, chemical):
    selected = unique if chemical is None else unique[unique.chemical == chemical]
    if selected.empty:
        raise ValueError("Requested chemical is absent from the unique-identity table")
    directory = OUT / "fresh_randomization"
    directory.mkdir(exist_ok=True)
    summary = []
    for row in selected.itertuples():
        targets = target_sets[row.chemical]
        bins = degree_bins(graph, 2 * max(len(targets), len(disease_set)))
        target_random = random_sets(graph, targets, bins, count, seed)
        disease_random = random_sets(graph, disease_set, bins, count, seed)
        null = []
        for iteration, (a, b) in enumerate(zip(target_random, disease_random), start=1):
            nearest = nx.multi_source_dijkstra_path_length(graph, b, weight=None)
            null.append({"iteration": iteration,
                         "closest": float(np.mean([nearest[n] for n in a])),
                         "n_targets": len(a), "n_disease": len(b)})
        frame = pd.DataFrame(null)
        slug = "".join(c for c in row.chemical if c.isalnum())
        frame.to_csv(directory / f"null_{slug}.csv.gz", index=False)
        mean = frame.closest.mean()
        std = frame.closest.std(ddof=0)
        zscore = (row.closest - mean) / std if std else math.nan
        empirical_p = (1 + (frame.closest <= row.closest).sum()) / (count + 1)
        summary.append({"chemical": row.chemical, "pubchem_cid": row.pubchem_cid,
                        "closest": row.closest, "mean_null": mean, "std_null": std,
                        "z_closest_new_run": zscore, "one_sided_empirical_p": empirical_p,
                        "n_random": count, "seed": seed})
        print(f"Randomization complete: {row.chemical} ({count} iterations)", flush=True)
    suffix = "all" if chemical is None else "".join(c for c in chemical if c.isalnum())
    pd.DataFrame(summary).to_csv(directory / f"summary_{suffix}.csv", index=False)
    metadata = {"seed": seed, "n_random": count, "python": platform.python_version(),
                "networkx": nx.__version__, "numpy": np.__version__,
                "disease_nodes": len(disease_set), "self_loops_retained": True,
                "sampling": "Upstream degree-bin sampling with sorted inputs; separate same-seed streams",
                "relation_to_historical_scores": "New run; historical node ordering and null vectors unavailable"}
    (directory / f"run_{suffix}.json").write_text(json.dumps(metadata, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["validate", "randomize"], default="validate")
    parser.add_argument("--n-random", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=452456)
    parser.add_argument("--chemical", help="Optional single chemical for a short smoke run")
    args = parser.parse_args()
    if args.n_random < 2:
        parser.error("--n-random must be at least 2")
    OUT.mkdir(parents=True, exist_ok=True)
    graph, disease, targets, saved, target_sets, disease_set = read_inputs()
    unique = original_records(targets, saved)
    check = recalculate_observed(graph, saved, target_sets, disease_set)
    ppi, centrality = ppi_topology()
    summary = {"source_proximity_records": len(saved), "unique_pubchem_identities": len(unique),
               "threshold_hit_records": int((saved.z_closest < -1.96).sum()),
               "threshold_hit_unique_identities": int((unique.z_closest < -1.96).sum()),
               "interactome_lcc_nodes": len(graph), "interactome_lcc_edges": graph.number_of_edges(),
               "interactome_self_loops": nx.number_of_selfloops(graph),
               "input_disease_genes": len(disease), "mapped_disease_genes": len(disease_set),
               "ppi_nodes": len(ppi), "ppi_edges": ppi.number_of_edges(),
               "ppi_component_sizes": sorted([len(c) for c in nx.connected_components(ppi)], reverse=True),
               "max_closest_absolute_error": float(check.closest_absolute_error.max()),
               "max_shortest_absolute_error": float(check.shortest_absolute_error.max()),
               "max_z_closest_error_from_saved_moments": float(check.z_closest_absolute_error.max()),
               "ppi_centrality_matches_original_cytoscape": True,
               "historical_z_scores": "Loaded from saved results; null vectors not present in source archive",
               "python": platform.python_version(), "networkx": nx.__version__,
               "numpy": np.__version__, "pandas": pd.__version__}
    (OUT / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.mode == "randomize":
        fresh_randomization(graph, unique, target_sets, disease_set,
                            args.n_random, args.seed, args.chemical)


if __name__ == "__main__":
    main()
