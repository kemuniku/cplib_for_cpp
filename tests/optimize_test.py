#!/usr/bin/env python3
"""最適化ビルドの実行、パス、探索順、失敗伝播と再コンパイルガード。"""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('optimize', ROOT / 'tools/optimize.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
compiler = sys.argv[1] if len(sys.argv) > 1 else 'g++'
with tempfile.TemporaryDirectory(prefix='cplib optimize ') as temporary:
    root = Path(temporary)
    first, second = root / 'first', root / 'second'
    first.mkdir(); second.mkdir()
    (first / 'chosen.hpp').write_text('#define CHOSEN 17\n')
    (second / 'chosen.hpp').write_text('#error wrong include priority\n')
    source, output = root / 'program with spaces.cpp', root / 'output with spaces'
    source.write_text('''#include <cplib/tmpl/optimize.hpp>
#include <chosen.hpp>
#if !defined(CPLIB_DEBUG)
static_assert(cplib::second_compile);
#else
static_assert(!cplib::second_compile);
#endif
int main(){return CHOSEN;}
''')
    assert module.optimize(source, output, compiler=compiler, include=[first,second]) == 0
    assert subprocess.run([str(output)], check=False).returncode == 17
    assert module.optimizeCpp(source, output, compiler=compiler, include=[first,second], debug=True) == 0
    assert module.optimize(source, output, compiler=compiler, include=[first,second], second_compile=True) == 0
    assert module.optimize(source, output, command=compiler+' -std=c++20 -DCPLIB_SECOND_COMPILE -O1', include=[first,second]) == 0
    try:
        module.optimize(source, output, command=compiler+' -O2')
        raise AssertionError('missing guard accepted')
    except ValueError:
        pass
    source.write_text('#error expected compilation failure\n')
    assert module.optimize(source, output, compiler=compiler) != 0
    qcfium = root / 'qcfium.cpp'
    qcfium.write_text('''#include <immintrin.h>
#include <cplib/tmpl/qcfium.hpp>
int main(){__m256i x=_mm256_set1_epi32(9);return _mm256_extract_epi32(x,0)-9;}
''')
    command=[compiler,'-std=c++20','-I',str(ROOT/'src'),str(qcfium),'-o',str(output)]
    if 'clang' in compiler:command += ['-mavx2','-O3','-funroll-loops']
    subprocess.run(command,check=True)
    subprocess.run([str(output)],check=True)
print('optimize tests passed')
