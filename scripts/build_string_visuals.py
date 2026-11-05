"""Prepare white-background STRING exports and verify their nodes and edges against the archived networks."""

import argparse
import hashlib
import io
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from PIL import Image, ImageChops, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "figures/reproduced/network/string_native"
SOURCES = ROOT / "data/network/string_visuals"
API = "https://version-12-0.string-db.org/api"


def svg_network(path):
    tree = ET.parse(path)
    nodes = {
        element.attrib["id"].removeprefix("node."): element.attrib["data-safe_div_label"]
        for element in tree.iter()
        if element.attrib.get("class") == "nwnodecontainer"
    }
    edges = set()
    for element in tree.iter():
        if element.attrib.get("class") != "nwlinkwrapper":
            continue
        _, first, second = element.attrib["id"].split(".")
        edges.add(tuple(sorted((nodes[first], nodes[second]))))
    return set(nodes.values()), edges


def white_background(source, destination):
    with Image.open(source) as image:
        image = image.convert("RGBA")
        white = Image.new("RGBA", image.size, "white")
        white.alpha_composite(image)
        white.convert("RGB").save(destination)


def notebook_preview(source, trim=False):
    destination = source.with_name(source.stem + "_preview.png")
    with Image.open(source) as image:
        image = image.convert("RGB")
        if trim:
            bounds = ImageChops.difference(image, Image.new("RGB", image.size, "white")).getbbox()
            if bounds:
                content = image.crop(bounds)
                image = Image.new("RGB", (content.width + 80, content.height + 80), "white")
                image.paste(content, (40, 40))
        image.thumbnail((2400, 2400), Image.Resampling.LANCZOS)
        image = image.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
        image.save(destination, optimize=True)
    return destination.relative_to(ROOT).as_posix()


def request(format_name, parameters):
    time.sleep(1.1)
    data = urllib.parse.urlencode(parameters).encode()
    query = urllib.request.Request(
        f"{API}/{format_name}/network", data=data,
        headers={"User-Agent": "CRC-research-archive/1.0"}
    )
    with urllib.request.urlopen(query, timeout=120) as response:
        payload = response.read()
    if format_name == "svg":
        ET.fromstring(payload)
    else:
        with Image.open(io.BytesIO(payload)) as image:
            image.verify()
    return payload


def panel_figure(entries, destination):
    font = ImageFont.load_default(size=38)
    small = ImageFont.load_default(size=27)
    panel_width, panel_height = 1500, 1390
    canvas = Image.new("RGB", (panel_width * 2, panel_height * 5), "white")
    draw = ImageDraw.Draw(canvas)
    for index, entry in enumerate(entries):
        x, y = index % 2 * panel_width, index // 2 * panel_height
        draw.text((x + 40, y + 22), f"{chr(65 + index)}. {entry['display_name']}", fill="#1d2529", font=font)
        draw.text((x + 40, y + 80),
                  f"{entry['nodes']} mapped proteins; {entry['edges']} associations", fill="#4b555a", font=small)
        with Image.open(ROOT / entry["image"]) as source:
            image = source.convert("RGB")
            image.thumbnail((panel_width - 65, panel_height - 150), Image.Resampling.LANCZOS)
            canvas.paste(image, (x + (panel_width - image.width) // 2,
                                 y + 135 + (panel_height - 150 - image.height) // 2))
    canvas.save(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--refresh", action="store_true")
    mode.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    SOURCES.mkdir(parents=True, exist_ok=True)
    nodes, edges = svg_network(ROOT / "figures/original/network/progression_ppi_string.svg")
    centrality = pd.read_csv(ROOT / "results/network/ppi_centrality_recalculated.csv")
    archived = pd.read_csv(ROOT / "data/network/string_interactions_short.tsv", sep="\t")
    archived_edges = {tuple(sorted(pair)) for pair in archived.iloc[:, :2].itertuples(index=False, name=None)}
    if nodes != set(centrality["name"]) or edges != archived_edges:
        raise ValueError("The archived STRING export differs from the analyzed PPI network")
    ppi_image = OUTPUT / "progression_ppi_string.png"
    white_background(ROOT / "figures/original/network/progression_ppi_string_high_resolution.png", ppi_image)
    ppi_preview = notebook_preview(ppi_image)
    with Image.open(ppi_image) as image:
        image.save(OUTPUT.parent / "figure_05A_ppi_network.png")
        image.save(OUTPUT.parent / "figure_05A_ppi_network.pdf", resolution=240)
    entries = []
    summary = ROOT / "results/network/compound_string_summary.csv"
    if not summary.exists():
        subprocess.run([sys.executable, str(ROOT / "scripts/python/fetch_string_compounds.py"),
                        "--offline"], cwd=ROOT, check=True)
    info = pd.read_csv(summary)
    for item in info.to_dict("records"):
        slug = "".join(letter for letter in item["chemical"] if letter.isalnum())
        mapping = pd.read_csv(ROOT / item["mapping_file"], sep="\t")
        identifiers = sorted(set(mapping["stringId"]))
        parameters = {
            "identifiers": "\r".join(identifiers), "species": 9606,
            "required_score": 400, "network_type": "functional",
            "network_flavor": "confidence", "add_nodes": 0,
            "add_color_nodes": 0, "add_white_nodes": 0,
            "hide_disconnected_nodes": 0, "hide_node_labels": 0,
            "block_structure_pics_in_bubbles": 0,
            "caller_identity": "crc-polyphenol-reproducibility"
        }
        metadata_path = SOURCES / f"{slug}.json"
        for format_name, suffix in (("svg", "svg"), ("highres_image", "png")):
            source = SOURCES / f"{slug}.{suffix}"
            if args.refresh or not source.exists():
                if args.offline:
                    parser.error(f"Missing cached STRING export: {source.relative_to(ROOT)}")
                source.write_bytes(request(format_name, parameters))
        labels, links = svg_network(SOURCES / f"{slug}.svg")
        expected_nodes = set(mapping["preferredName"])
        table = pd.read_csv(ROOT / item["edges_file"], sep="\t")
        expected_edges = {tuple(sorted(pair)) for pair in
                          table.loc[table["score"] >= 0.4,
                                    ["preferredName_A", "preferredName_B"]].itertuples(index=False, name=None)}
        if labels != expected_nodes or links != expected_edges:
            raise ValueError(f"STRING image for {item['chemical']} differs from its archived edge table: "
                             f"{len(labels)}/{len(expected_nodes)} nodes; {len(links)}/{len(expected_edges)} edges")
        white_background(SOURCES / f"{slug}.png", OUTPUT / f"{slug}.png")
        preview = notebook_preview(OUTPUT / f"{slug}.png", trim=len(labels) <= 4)
        if args.refresh or not metadata_path.exists():
            metadata = {"api_base": API, "parameters": parameters,
                        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                        "svg_sha256": hashlib.sha256((SOURCES / f"{slug}.svg").read_bytes()).hexdigest(),
                        "png_sha256": hashlib.sha256((SOURCES / f"{slug}.png").read_bytes()).hexdigest()}
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        entry = {"chemical": item["chemical"], "display_name": item["display_name"],
                 "nodes": len(labels), "edges": len(links),
                 "image": (OUTPUT / f"{slug}.png").relative_to(ROOT).as_posix(),
                 "notebook_image": preview,
                 "source_svg": (SOURCES / f"{slug}.svg").relative_to(ROOT).as_posix(),
                 "source_png": (SOURCES / f"{slug}.png").relative_to(ROOT).as_posix(),
                 "request_metadata": metadata_path.relative_to(ROOT).as_posix(),
                 "nodes_match_archive": True, "edges_match_archive": True}
        entries.append(entry)
        print(f"{item['display_name']}: {len(labels)} nodes; {len(links)} edges", flush=True)
    panel_figure(entries, OUTPUT / "compound_target_networks_string.png")
    with Image.open(OUTPUT / "compound_target_networks_string.png") as image:
        image.save(OUTPUT.parent / "figure_08_compound_target_networks_reconstructed.png")
        image.save(OUTPUT.parent / "figure_08_compound_target_networks_reconstructed.pdf", resolution=240)
    manifest = {"preview_processing": {"maximum_long_edge_px": 2400,
                                       "resampling": "Lanczos",
                                       "palette_colors": 256,
                                       "background": "white",
                                       "small_network_trimming": "For networks of up to four nodes, crop to all nonwhite pixels and add 40 pixels of white padding on every side",
                                       "full_resolution_rgb_preserved": True},
                "ppi": {"image": ppi_image.relative_to(ROOT).as_posix(),
                        "notebook_image": ppi_preview,
                        "source": "figures/original/network/progression_ppi_string_high_resolution.png",
                        "nodes": len(nodes), "edges": len(edges),
                        "nodes_match_archive": True, "edges_match_archive": True},
                "compounds": entries}
    (ROOT / "results/network/string_visuals.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
