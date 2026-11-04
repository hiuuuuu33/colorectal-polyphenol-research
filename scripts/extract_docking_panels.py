"""Extract the four original docking panels from the archived PowerPoint without altering their image bytes."""

from pathlib import Path
from zipfile import ZipFile
import csv
import hashlib
import io
import posixpath
import xml.etree.ElementTree as ET

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
PRESENTATION = ROOT / "data/docking/records/docking_results.pptx"
OUTPUT = ROOT / "figures/original/docking"
NAMESPACES = {
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "p": "http://schemas.openxmlformats.org/presentationml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
PANELS = (
    (1, "CXCL8-phloretin", ("phloretin_6wzm_3d.png", "phloretin_6wzm_interactions.png")),
    (2, "ADH1C-kaempferol", ("kaempferol_1ht0_3d.png", "kaempferol_1ht0_interactions.png")),
)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = []
    with ZipFile(PRESENTATION) as archive:
        for slide_number, expected_title, filenames in PANELS:
            slide_path = f"ppt/slides/slide{slide_number}.xml"
            slide = ET.fromstring(archive.read(slide_path))
            title = "".join(slide.itertext()).strip()
            if title != expected_title:
                raise ValueError(f"Unexpected title on slide {slide_number}: {title}")
            relationships = ET.fromstring(archive.read(
                f"ppt/slides/_rels/slide{slide_number}.xml.rels"
            ))
            targets = {item.attrib["Id"]: item.attrib["Target"] for item in relationships}
            pictures = slide.findall(".//p:pic/p:blipFill/a:blip", NAMESPACES)
            if len(pictures) != len(filenames):
                raise ValueError(f"Expected two original images on slide {slide_number}")
            for picture, filename in zip(pictures, filenames):
                relationship_id = picture.attrib[f"{{{NAMESPACES['r']}}}embed"]
                media = posixpath.normpath(posixpath.join(
                    posixpath.dirname(slide_path), targets[relationship_id]
                ))
                image_bytes = archive.read(media)
                with Image.open(io.BytesIO(image_bytes)) as image:
                    if image.format != "PNG":
                        raise ValueError(f"Expected an original PNG image: {media}")
                    width, height = image.size
                    image.verify()
                destination = OUTPUT / filename
                destination.write_bytes(image_bytes)
                rows.append({
                    "file": destination.relative_to(ROOT).as_posix(),
                    "source_pptx": PRESENTATION.relative_to(ROOT).as_posix(),
                    "slide": slide_number,
                    "slide_title": title,
                    "ppt_media": media,
                    "width": width,
                    "height": height,
                    "sha256": hashlib.sha256(image_bytes).hexdigest(),
                })
    manifest = ROOT / "metadata/docking_panels.csv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Extracted {len(rows)} original docking panels without changing their image bytes.")


if __name__ == "__main__":
    main()
