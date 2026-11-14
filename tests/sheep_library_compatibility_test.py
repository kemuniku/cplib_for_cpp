"""Check library headers with the unchanged sheep template in both include orders."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BODY = '''
using namespace cplib;
#ifndef len
#error sheep len macro must remain available
#endif
#ifndef read
#error sheep read macro must remain available
#endif
#ifndef out
#error sheep out macro must remain available
#endif
int main() {
    vector<int> v; set<int> s; map<int,int> m; hash<int> h;
    queue<int> q; stack<int> st; list<int> l; deque<int> d;
    array<int,1> a; pair<int,int> p; tuple<int> t; priority_queue<int> pq;
    bitset<8> b; optional<int> o; variant<int> va; function<int()> f;
    string text; unordered_set<int> us; unordered_map<int,int> um; complex<double> c;
    return len(v);
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    headers = sorted(p.relative_to(ROOT / 'src').as_posix()
                     for p in (ROOT / 'src/cplib').rglob('*.hpp')
                     if 'tmpl' not in p.relative_to(ROOT / 'src/cplib').parts)
    cases = [(header, order) for header in headers for order in (0, 1)]

    def check(case):
        header, order = case
        includes = ['sheep_library_macros.hpp', header]
        if order:
            includes.reverse()
        source = ''.join(f'#include <{h}>\n' for h in includes) + BODY
        p = subprocess.run([args.compiler, '-std=gnu++20', '-fsyntax-only',
                            '-fmax-errors=5', '-I', str(ROOT / 'src'), '-I', str(ROOT / 'tests'), '-x', 'c++', '-'],
                           input=source, text=True, capture_output=True, timeout=60)
        return case, p.returncode, p.stderr

    failed = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = [pool.submit(check, case) for case in cases]
        for i, future in enumerate(as_completed(pending), 1):
            case, code, errors = future.result()
            if code:
                failed += 1
                print(f'FAIL: {case}\n{errors}', flush=True)
            if i % 100 == 0 or i == len(cases):
                print(f'{i}/{len(cases)} checked; {failed} failed', flush=True)
    raise SystemExit(bool(failed))


if __name__ == '__main__':
    main()
