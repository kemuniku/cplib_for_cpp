#!/usr/bin/env python3
"""各ヘッダーを単独includeし、2翻訳単位でリンクできることを検証する。

別実装は同じ cplib の公開名を持つため、ヘッダーごとに独立した実行ファイルを作る。
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    root = Path(__file__).resolve().parents[1]
    headers = sorted((root / 'src/cplib').rglob('*.hpp'))
    command = [args.compiler, '-std=c++20', '-O0', '-I', str(root / 'src')]

    def check(header):
        name = header.relative_to(root / 'src')
        with tempfile.TemporaryDirectory(prefix='cppcplib_headers_') as directory:
            temp = Path(directory)
            first, second = temp / 'first.cpp', temp / 'second.cpp'
            include = f'#include <{name}>\n'
            # Each translation unit must compile independently, and the same
            # header must not produce duplicate symbols when linked.
            first.write_text(include + 'int helper(){return 0;}\n')
            second.write_text(include + 'int helper(); int main(){return helper();}\n')
            result = subprocess.run(command + [str(first), str(second), '-o', str(temp / 'check')],
                                    capture_output=True, text=True)
            return name, result

    failures = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for index, (name, result) in enumerate(pool.map(check, headers), 1):
            if result.returncode:
                failures += 1
                print(f'{name}:\n{result.stderr}', flush=True)
            if index % 50 == 0:
                print(f'{index}/{len(headers)} headers checked', flush=True)
    if failures:
        raise SystemExit(f'{failures} headers failed')
    print(f'{len(headers)} headers: standalone and multi-TU checks passed')


if __name__ == '__main__':
    main()
