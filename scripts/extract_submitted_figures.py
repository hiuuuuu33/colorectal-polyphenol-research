"""Extract figures from a supplied local PDF, recording checksums and page numbers."""
from pathlib import Path
import argparse
import csv
import hashlib
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures/submitted'
PAGES = {1: 8, 2: 9, 3: 10, 4: 12, 5: 13, 6: 15, 7: 18, 8: 19, 9: 21}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', type=Path, required=True, help='Local manuscript PDF')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    doc = pymupdf.open(args.pdf)
    rows = []
    for number, page_number in PAGES.items():
        page = doc[page_number - 1]
        images = page.get_images(full=True)
        if len(images) != 1:
            raise ValueError(f'Expected one embedded figure on page {page_number}')
        extracted = doc.extract_image(images[0][0])
        path = OUT / f'figure_{number:02d}.{extracted["ext"]}'
        path.write_bytes(extracted['image'])
        rows.append({'figure': number, 'pdf_page': page_number,
                     'file': path.relative_to(ROOT).as_posix(),
                     'width': extracted['width'], 'height': extracted['height'],
                     'sha256': hashlib.sha256(extracted['image']).hexdigest()})
    with (ROOT / 'metadata/submitted_figures.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'Extracted {len(rows)} embedded figures without changing the PDF.')


if __name__ == '__main__':
    main()
