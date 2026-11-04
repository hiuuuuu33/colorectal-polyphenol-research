"""Build Figure 7 from cached PubChem structures for the top ten compounds."""
from pathlib import Path
import argparse
import csv
import datetime
import io
import json
import subprocess
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true', help='Fetch missing PubChem PNGs.')
    args = parser.parse_args()
    cache = ROOT / 'data/chemistry/pubchem'
    cache.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((ROOT / 'results/network/top10_candidates.csv').open()))
    canvas = Image.new('RGB', (1400, 2250), 'white')
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype('DejaVuSans.ttf', 25)
        small = ImageFont.truetype('DejaVuSans.ttf', 21)
    except OSError:
        font = ImageFont.load_default(size=25)
        small = ImageFont.load_default(size=21)
    for i, row in enumerate(rows):
        cid = row['pubchem_cid']
        filename = cache / f'CID_{cid}.png'
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size=large'
        if not filename.exists():
            if not args.download:
                raise FileNotFoundError(f'{filename.name}: run once with --download')
            content = subprocess.run(
                ['curl', '--fail', '--location', '--silent', '--show-error',
                 '--retry', '2', '--connect-timeout', '15', '--max-time', '60', url],
                check=True, capture_output=True).stdout
            Image.open(io.BytesIO(content)).verify()
            filename.write_bytes(content)
            (cache / f'CID_{cid}.json').write_text(json.dumps({
                'cid': cid, 'source_url': url,
                'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()
            }, indent=2))
        image = Image.open(filename).convert('RGB')
        image.thumbnail((620, 330))
        x, y = (i % 2) * 700, (i // 2) * 450
        canvas.paste(image, (x + (700-image.width)//2, y+60+(330-image.height)//2))
        name = row['display_name'].replace('naphtoquinone', 'naphthoquinone')
        draw.text((x+28, y+12), f'{chr(65+i)}  {name}', font=font, fill='#202020')
        draw.text((x+28, y+405), f'PubChem CID {cid}', font=small, fill='#555555')
    output = ROOT / 'figures/reproduced/chemistry/figure_07_compound_structures.png'
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output)
    print(f'Built {output.relative_to(ROOT)} from 10 cached PubChem depictions.')


if __name__ == '__main__':
    main()
