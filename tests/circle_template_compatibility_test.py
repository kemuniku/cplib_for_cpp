"""円・最小包含円と既存テンプレートの同時使用を、include順も変えて確認する。"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import permutations
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HEADERS = ['cplib/geometry/circle.hpp', 'cplib/geometry/circle_integer_impl.hpp',
           'cplib/geometry/minimum_enclosing_circle.hpp', 'cplib/str/lcs.hpp',
           'cplib/math/int128.hpp', 'cplib/math/fractions.hpp']
BODY = '''
int main() {
    ll size = 3;
    vector<Int> values{1, 2, 3};
    rep(i, size) { assert(values[i] == i + 1); }
    auto c = initCircle(Point<Int>{0, 0}, Int(2));
    auto m = minimum_enclosing_circle(vector<Point<Int>>{{0, 0}, {1, 0}});
    assert(c.contains(Point<Int>{1, 0}));
    assert(c.intersection_count(Line<Int>{{-3, 0}, {3, 0}}) == 2);
    assert(c.cross_points_approx(Line<Int>{{-3, 0}, {3, 0}}).points.size() == 2);
    assert(c.tangent_lines_approx(Point<Int>{3, 0}).tangents.size() == 2);
    assert(m.radius_squared_exact() == Fraction<Int>(1, 4));
    assert(LCS(values, vector<Int>{2, 3}) == 2);
    assert(cplib::to_string(parseInt128("-170141183460469231731687303715884105728"))
           == "-170141183460469231731687303715884105728");
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--syntax-only', action='store_true', help='実行せず構文検査だけを行う')
    args = parser.parse_args()
    templates = ['cplib/tmpl/sheep.hpp', 'cplib/tmpl/compact.hpp']
    cases = []
    for order in permutations(templates):
        for library_first in (False, True):
            for using_first in (False, True):
                cases.append((order, library_first, using_first))
    with tempfile.TemporaryDirectory(prefix='circle_templates_') as directory:
        def check(index_case):
            index, (order, library_first, using_first) = index_case
            prefix = HEADERS if library_first else list(order)
            suffix = list(order) if library_first else HEADERS
            source = ''.join(f'#include <{header}>\n' for header in prefix)
            if using_first:
                source += 'using namespace cplib;\n'
            source += ''.join(f'#include <{header}>\n' for header in suffix)
            if not using_first:
                source += 'using namespace cplib;\n'
            source += BODY
            file = Path(directory) / f'case_{index}.cpp'
            binary = Path(directory) / f'case_{index}'
            file.write_text(source)
            command = [args.compiler, '-std=c++20', '-O1', '-I', str(ROOT / 'src'), str(file)]
            command += ['-fsyntax-only'] if args.syntax_only else ['-o', str(binary)]
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise AssertionError(f'{order}, library_first={library_first}, using_first={using_first}\n{result.stderr}')
            if not args.syntax_only:
                subprocess.run([str(binary)], check=True)
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(check, enumerate(cases)))
    mode = 'compile' if args.syntax_only else 'runtime'
    print(f'{len(cases)} template/include-order {mode} checks passed')


if __name__ == '__main__':
    main()
