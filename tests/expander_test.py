#!/usr/bin/env python3
"""提出ファイル単体でのコンパイル、重複・循環・条件付きincludeを検証する。"""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('expander', ROOT / 'expander.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
compiler = sys.argv[1] if len(sys.argv) > 1 else 'g++'


def run_expander(source, *args, cwd):
    return subprocess.run([sys.executable, str(ROOT / 'expander.py'), str(source), *args],
                          cwd=cwd, text=True, capture_output=True, check=False)


def compile_and_run(source, output, *flags):
    # -I を指定せず、生成ファイルだけで提出できることを検証する。
    subprocess.run([compiler, '-std=c++20', *flags, str(source), '-o', str(output)], check=True)
    subprocess.run([str(output)], check=True)


with tempfile.TemporaryDirectory(prefix='cplib expander ') as temporary:
    root = Path(temporary)
    source = root / 'main with spaces.cpp'
    output = root / 'submission.cpp'
    binary = root / 'submission'
    source.write_text('''#include <cplib/math/floor_sum.hpp>
#include "cplib/collections/unionfind.hpp"
#include <cplib/math/../math/floor_sum.hpp>
#include <cassert>
int main() {
    assert(cplib::floor_sum(4, 10, 6, 3) == 3);
    cplib::UnionFind uf(3);
    uf.unite(0, 2);
    assert(uf.issame(0, 2) && !uf.issame(0, 1));
}
''', encoding='utf-8')
    original = source.read_text()
    result = run_expander(source, cwd=root)
    assert result.returncode == 0, result.stderr
    assert not result.stderr
    expanded = result.stdout
    assert '#include <cplib/' not in expanded
    assert '#include "cplib/' not in expanded
    assert '#pragma once' not in expanded
    assert expanded.count('using Int = std::int64_t;') == 1
    assert expanded.count('inline Int floor_sum(') == 1
    assert '#include <cassert>' in expanded
    result = run_expander(source, '-o', str(output), cwd=root)
    assert result.returncode == 0, result.stderr
    assert not result.stdout
    assert output.read_text() == expanded
    assert source.read_text() == original
    compile_and_run(output, binary)
    assert run_expander(source, '-o', '-', cwd=root).stdout == expanded

    # 同一ヘッダへの別表記、symlink、ダイヤ形依存、循環参照。
    first, second = root / 'first', root / 'second'
    first.mkdir()
    second.mkdir()
    common = first / 'common.hpp'
    common.write_text('#pragma once\nstruct Shared {};\n', encoding='utf-8')
    (first / 'alias.hpp').symlink_to(common)
    (first / 'left.hpp').write_text('#pragma once\n#include "common.hpp"\n#include "right.hpp"\n', encoding='utf-8')
    (first / 'right.hpp').write_text('#pragma once\n#include "./common.hpp"\n#include "left.hpp"\n', encoding='utf-8')
    (second / 'common.hpp').write_text('#error wrong include search order\n', encoding='utf-8')
    source.write_text('#include <left.hpp>\n#include <right.hpp>\n#include <alias.hpp>\nint main() { Shared s; }\n', encoding='utf-8')
    result = run_expander(source, '-I', str(first), '-I', str(second), cwd=root)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count('struct Shared {};') == 1
    output.write_text(result.stdout, encoding='utf-8')
    compile_and_run(output, binary)

    # 最初のincludeが無効な分岐でも、後続の定義を失わない。
    source.write_text('''#if 0
#include "first/common.hpp"
#endif
#if defined(SELECT_FIRST)
#include "first/left.hpp"
#elif defined(SELECT_SECOND)
#include "first/right.hpp"
#else
#include "first/alias.hpp"
#endif
#include "first/common.hpp"
int main() { Shared s; }
''', encoding='utf-8')
    output.write_text(module.Expander().expand(source), encoding='utf-8')
    for flags in [(), ('-DSELECT_FIRST',), ('-DSELECT_SECOND',)]:
        compile_and_run(output, binary, *flags)

    # コメント・raw stringのincludeを無視し、実際のinclude行のコメントは保持。
    source.write_text('''/*
#include "missing-in-comment.hpp"
*/
// continued comment \\
#include "missing-in-continued-comment.hpp"
const char* text = R"tag(
#include "missing-in-string.hpp"
)tag";
/* leading
*/ # include \\
"first/common.hpp" /* trailing
comment */
int main() { Shared s; }
''', encoding='utf-8')
    output.write_text(module.Expander().expand(source), encoding='utf-8')
    assert 'missing-in-comment.hpp' in output.read_text()
    assert 'missing-in-string.hpp' in output.read_text()
    compile_and_run(output, binary)

    # 見つからないローカルヘッダは失敗し、既存の出力や入力を壊さない。
    for include in ['<cplib/does_not_exist.hpp>', '"does_not_exist.hpp"']:
        source.write_text(f'#include {include}\n', encoding='utf-8')
        output.write_text('keep existing output\n', encoding='utf-8')
        result = run_expander(source, '-o', str(output), cwd=root)
        assert result.returncode != 0
        assert 'does_not_exist.hpp' in result.stderr
        assert not result.stdout
        assert output.read_text() == 'keep existing output\n'
    source.write_text('#include "first/common.hpp"\n', encoding='utf-8')
    for destination in [source, common]:
        before = destination.read_text()
        result = run_expander(source, '-o', str(destination), cwd=root)
        assert result.returncode != 0
        assert destination.read_text() == before
    result = run_expander(root / 'missing.cpp', cwd=root)
    assert result.returncode != 0
    assert not result.stdout

print('expander tests passed')
