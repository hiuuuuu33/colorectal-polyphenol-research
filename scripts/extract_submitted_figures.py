"""Extract embedded figures from the unchanged submitted PDF."""
from pathlib import Path
import csv
import hashlib
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'manuscript/submitted-manuscript.pdf'
OUT = ROOT / 'figures/submitted'
PAGES = {1: 8, 2: 9, 3: 10, 4: 12, 5: 13, 6: 15, 7: 18, 8: 19, 9: 21}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    doc = pymupdf.open(PAPER)
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
    (ROOT / 'manuscript/submitted-manuscript.txt').write_text(
        '\n\n'.join(f'PDF PAGE {i + 1}\n{p.get_text()}' for i, p in enumerate(doc)),
        encoding='utf-8')
    print(f'Extracted {len(rows)} embedded figures without changing the PDF.')


if __name__ == '__main__':
    main()
