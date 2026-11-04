"""Verify SHA-256 checksums for research inputs and optionally update the manifest."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'metadata/input_checksums.csv'
IMMUTABLE = ('data', 'archive', 'upstream', 'manuscript', 'figures/submitted', 'figures/original')


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-manifest', action='store_true')
    args = parser.parse_args()
    if args.write_manifest:
        rows = []
        for name in IMMUTABLE:
            for path in sorted((ROOT / name).rglob('*')):
                if path.is_file() and not any(p.startswith('.') for p in path.relative_to(ROOT).parts):
                    rows.append({'path': path.relative_to(ROOT).as_posix(),
                                 'bytes': path.stat().st_size, 'sha256': sha(path)})
        with MANIFEST.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['path', 'bytes', 'sha256'])
            writer.writeheader()
            writer.writerows(rows)
    rows = list(csv.DictReader(MANIFEST.open()))
    failed = [row['path'] for row in rows if not (ROOT / row['path']).is_file()
              or sha(ROOT / row['path']) != row['sha256']]
    print(json.dumps({'checked_input_files': len(rows), 'checksum_failures': failed}, indent=2))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
