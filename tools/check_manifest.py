#!/usr/bin/env python3
"""移植対象の抜け・元ソース変更・対応ファイルの欠落を検査する。"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1] / 'reference/cplib')
parser.add_argument('--require-complete', action='store_true')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
entries = json.loads((root / 'port_manifest.json').read_text())
listed = {e['source'] for e in entries}
actual = {str(p.relative_to(args.source)) for p in args.source.rglob('*') if p.is_file()}
if args.source.resolve() == (root / 'reference/cplib').resolve():
    actual |= {e['source'] for e in entries if e['status'] == 'preserved_binary'}
assert listed == actual, f'missing={actual-listed}, extra={listed-actual}'
assert len(listed) == len(entries), 'duplicate entries'
for entry in entries:
    source = args.source / entry['source']
    if entry['status'] == 'preserved_binary' and not source.exists():
        source = root / entry['target']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == entry['sha256'], f'source changed: {source}'
    if entry['status'] != 'pending':
        assert (root / entry['target']).is_file(), f'missing target: {entry["target"]}'
    if entry['status'] == 'preserved_binary':
        assert (root / entry['target']).read_bytes() == source.read_bytes()
counts = Counter(e['status'] for e in entries)
print(dict(counts))
if args.require_complete:
    assert all(e['status'] in ('ported', 'preserved_binary') for e in entries), 'port incomplete'
