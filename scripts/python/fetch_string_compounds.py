#!/usr/bin/env python3
"""Build compound target networks from STRING 12.0 mappings and edges for Figure 8.

Use --offline to rebuild summaries and layouts from cached responses, or
--refresh to retrieve new responses from STRING.
"""

import argparse
import csv
import hashlib
import io
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "network" / "string_compounds"
API = "https://version-12-0.string-db.org/api/tsv"


def readable_layout(graph):
    if len(graph) < 10:
        return nx.circular_layout(graph)
    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    positions = dict(nx.spring_layout(graph.subgraph(components[0]),
                                     seed=2026, iterations=250, scale=1))
    for index, component in enumerate(components[1:]):
        center_x = -1 + (index % 9) * 0.25
        center_y = -1.35 - (index // 9) * 0.26
        layout = nx.circular_layout(graph.subgraph(component), scale=0.065)
        for node, (x, y) in layout.items():
            positions[node] = (center_x + x, center_y + y)
    return positions


def request(endpoint, parameters):
    time.sleep(1.05)
    parameters = dict(parameters, species=9606, caller_identity="crc-polyphenol-reproducibility")
    body = urllib.parse.urlencode(parameters).encode()
    req = urllib.request.Request(f"{API}/{endpoint}", data=body,
                                 headers={"User-Agent": "CRC-research-archive/1.0"})
    with urllib.request.urlopen(req, timeout=90) as response:
        payload = response.read().decode("utf-8")
    if payload.lstrip().startswith("<"):
        raise RuntimeError("STRING returned an HTML response")
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--refresh", action="store_true")
    mode.add_argument("--offline", action="store_true",
                      help="Rebuild summaries and layouts using only existing cache files")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(ROOT / "results/network/top10_recorded_targets.csv")
    if args.offline:
        missing = []
        for name in targets.chemical.drop_duplicates():
            slug = "".join(character for character in name if character.isalnum())
            for suffix in ("mapping.tsv", "edges.tsv", "request.json"):
                path = OUT / f"{slug}_{suffix}"
                if not path.is_file():
                    missing.append(path.relative_to(ROOT).as_posix())
        if missing:
            parser.error("Offline cache files are missing:\n" + "\n".join(missing))
    entries, all_positions, manifests = [], [], []
    for name, subset in targets.groupby("chemical", sort=False):
        slug = "".join(character for character in name if character.isalnum())
        mapping_file = OUT / f"{slug}_mapping.tsv"
        edge_file = OUT / f"{slug}_edges.tsv"
        metadata_file = OUT / f"{slug}_request.json"
        identifiers = sorted(set(subset.symbol.dropna().astype(str)))
        if args.refresh or not mapping_file.exists():
            if args.offline:
                parser.error(f"Offline cache file is missing: {mapping_file.relative_to(ROOT)}")
            payload = request("get_string_ids", {"identifiers": "\r".join(identifiers),
                                                 "limit": 1, "echo_query": 1})
            mapping_file.write_text(payload)
        mapping = pd.read_csv(mapping_file, sep="\t")
        if "stringId" not in mapping.columns:
            raise RuntimeError(f"Unresolved STRING mapping for {name}: {mapping_file}")
        string_ids = sorted(set(mapping.stringId))
        if args.refresh or not edge_file.exists():
            if args.offline:
                parser.error(f"Offline cache file is missing: {edge_file.relative_to(ROOT)}")
            if len(string_ids) > 1:
                payload = request("network", {"identifiers": "\r".join(string_ids),
                                              "required_score": 400, "network_type": "functional",
                                              "add_nodes": 0})
                if not payload.strip():
                    payload = "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\tnscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n"
            else:
                payload = "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\tnscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n"
            edge_file.write_text(payload)
        edges = pd.read_csv(edge_file, sep="\t")
        graph = nx.Graph()
        graph.add_nodes_from(string_ids)
        if len(edges):
            if "score" not in edges:
                raise RuntimeError(f"Unexpected STRING edge response for {name}")
            edges = edges[(edges.stringId_A.isin(string_ids)) &
                          (edges.stringId_B.isin(string_ids)) & (edges.score >= 0.4)]
            graph.add_edges_from(zip(edges.stringId_A, edges.stringId_B))
        label_lookup = dict(zip(mapping.stringId, mapping.preferredName))
        pos = readable_layout(graph)
        for string_id in graph:
            x, y = pos[string_id]
            all_positions.append({"chemical": name, "rank": int(subset["rank"].iloc[0]),
                                  "string_id": string_id, "symbol": label_lookup[string_id],
                                  "x": x, "y": y, "degree": graph.degree(string_id)})
        entries.append({"rank": int(subset["rank"].iloc[0]), "chemical": name,
                        "display_name": subset.display_name.iloc[0],
                        "recorded_targets": len(subset), "string_mapped_nodes": len(graph),
                        "edges_score_ge_0_4": graph.number_of_edges(),
                        "mapping_file": mapping_file.relative_to(ROOT).as_posix(),
                        "edges_file": edge_file.relative_to(ROOT).as_posix()})
        if args.refresh or not metadata_file.exists():
            if args.offline:
                parser.error(f"Offline cache file is missing: {metadata_file.relative_to(ROOT)}")
            metadata = {"chemical": name, "api_base": API, "string_version": "12.0",
                        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                        "species": 9606, "network_type": "functional", "required_score": 400,
                        "add_nodes": 0, "singleton_rule": "Keep single query node without added neighbors",
                        "query_symbols": identifiers,
                        "mapping_sha256": hashlib.sha256(mapping_file.read_bytes()).hexdigest(),
                        "edges_sha256": hashlib.sha256(edge_file.read_bytes()).hexdigest()}
            metadata_file.write_text(json.dumps(metadata, indent=2) + "\n")
        elif args.offline:
            json.loads(metadata_file.read_text())
        print(f"{name}: {len(subset)} targets; {len(graph)} mapped nodes; {graph.number_of_edges()} edges", flush=True)
    pd.DataFrame(entries).to_csv(ROOT / "results/network/compound_string_summary.csv", index=False)
    pd.DataFrame(all_positions).to_csv(ROOT / "results/network/compound_string_layouts.csv", index=False)


if __name__ == "__main__":
    main()
