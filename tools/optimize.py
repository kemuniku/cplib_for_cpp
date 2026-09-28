#!/usr/bin/env python3
"""tmpl/optimize.nim の C++ ビルド側移植。シェルを介さずコンパイラを実行する。"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MARKER = '-DCPLIB_SECOND_COMPILE'


def optimize(source, output, *, compiler=None, command=None, include=(), flags=(), debug=False, second_compile=False):
    """最適化コンパイル、出力先・探索順の保持、失敗コードの伝播を行う。"""
    source = Path(source).resolve()
    output = Path(output).resolve()
    includes = [str(Path(p).resolve()) for p in include]
    if command and not debug and not second_compile:
        argv = shlex.split(command)
        if not any(a == MARKER or a.startswith(MARKER + '=') for a in argv):
            raise ValueError('custom command must define CPLIB_SECOND_COMPILE')
    else:
        argv = shlex.split(compiler or os.environ.get('CXX', 'g++'))
        argv += ['-std=c++20']
        if debug:
            argv += ['-DCPLIB_DEBUG', '-g']
        elif not second_compile:
            argv += [MARKER, '-DNDEBUG', '-O3', '-flto', '-m64', '-march=native',
                     '-ffast-math', '-funroll-loops']
            # ClangにはGCC固有のIPA-PTAが存在しない。
            version = subprocess.run(argv[:1] + ['--version'], capture_output=True, text=True, check=False)
            if 'clang' not in (version.stdout + version.stderr).lower():
                argv += ['-fipa-pta']
        else:
            argv += [MARKER]
    argv += ['-I', str(ROOT / 'src')]
    for path in includes:
        argv += ['-I', path]
    argv += list(flags) + ['-o', str(output), str(source)]
    print('--- Compiling with ' + ('debug' if debug else 'optimized') + ' settings ---', flush=True)
    print(shlex.join(argv), flush=True)
    return subprocess.run(argv, cwd=source.parent, check=False).returncode


# NimのC/C++バックエンドの区別は、移植後はいずれもC++コンパイラに対応する。
optimizeCpp = optimize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('-o', '--output', default='a.out')
    parser.add_argument('--compiler')
    parser.add_argument('--command', help='custom compiler command; must contain -DCPLIB_SECOND_COMPILE')
    parser.add_argument('-I', '--include', action='append', default=[])
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--second-compile', action='store_true')
    arguments = sys.argv[1:]
    extra = []
    if '--' in arguments:
        i = arguments.index('--')
        arguments, extra = arguments[:i], arguments[i + 1:]
    args = parser.parse_args(arguments)
    try:
        return optimize(args.source, args.output, compiler=args.compiler, command=args.command,
                        include=args.include, flags=extra, debug=args.debug, second_compile=args.second_compile)
    except (ValueError, OSError) as error:
        print(error, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
