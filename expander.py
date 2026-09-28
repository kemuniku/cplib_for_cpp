#!/usr/bin/env python3
"""cppcplib の include を展開し、ジャッジ提出用の C++ ソースを生成する。"""

import argparse
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parent
# コメント・raw string 内の #include は処理しない。位置と改行数は維持する。
TOKENS = re.compile(
    r'//(?:\\\n|[^\n])*|/\*[\s\S]*?\*/'
    r'|R"(?P<delimiter>[^ ()\\\t\r\n]{0,16})\([\s\S]*?\)(?P=delimiter)"'
    r'|"(?:\\[\s\S]|[^"\\])*"'
    r"|'(?:\\[^\n]|[^'\\\n])*'"
)
DIRECTIVE = re.compile(
    r'^[ \t\v\f\r]*#[ \t]*(?P<name>\w+)\b(?P<rest>(?:\\\n|[^\n])*)',
    re.MULTILINE,
)
INCLUDE = re.compile(r'(?:\s|\\\n)*(?:<(?P<angle>[^>\n]+)>|"(?P<quote>[^"\n]+)")')
ONCE = re.compile(r'\s+once\b')


def directive_text(source):
    """解析用のコピーだけからコメントと複数行文字列を隠す。"""
    def mask(match):
        token = match.group()
        if token.startswith(('//', '/*', 'R"')) or '\n' in token:
            return re.sub(r'[^\n]', ' ', token)
        return token

    return TOKENS.sub(mask, source)


class Expander:
    def __init__(self, include=()):
        self.include_dirs = [Path(p).resolve() for p in include] + [ROOT / 'src']
        self.files = {}

    def resolve_include(self, name, quoted, parent):
        directories = ([parent] if quoted else []) + self.include_dirs
        for directory in directories:
            candidate = directory / name
            if candidate.is_file():
                return candidate.resolve()
        if quoted or name.startswith('cplib/'):
            raise ValueError(f'include が見つかりません: {name}')
        # 標準ヘッダやジャッジ側の外部ライブラリはそのまま残す。
        return None

    def expand(self, source):
        self.files.clear()
        return self._expand(Path(source).resolve(), set(), root=True)

    def _expand(self, path, seen, *, root=False):
        # resolve 済みパスで比較し、別表記・symlink・循環参照も重複させない。
        if path in seen:
            return ''
        seen.add(path)
        self.files.setdefault(path, len(self.files))
        source = path.read_text(encoding='utf-8-sig')
        if source and not source.endswith('\n'):
            source += '\n'
        masked = directive_text(source)
        output = []
        cursor = 0
        branches = []
        for directive in DIRECTIVE.finditer(masked):
            kind = directive['name']
            rest = directive['rest']
            if kind in ('if', 'ifdef', 'ifndef'):
                branches.append((seen.copy(), [], False))
            elif kind in ('elif', 'else') and branches:
                before, completed, has_else = branches[-1]
                completed.append(seen.copy())
                seen.clear()
                seen.update(before)
                branches[-1] = (before, completed, has_else or kind == 'else')
            elif kind == 'endif' and branches:
                before, completed, has_else = branches.pop()
                completed.append(seen.copy())
                if not has_else:
                    completed.append(before)
                # 全分岐で読み込まれたファイルだけを展開済みとする。
                seen.intersection_update(*completed)

            replacement = None
            end = directive.end()
            if kind == 'include':
                match = INCLUDE.match(rest)
                if match:
                    name = match['angle'] or match['quote']
                    try:
                        child = self.resolve_include(name, match['quote'] is not None, path.parent)
                        if child is not None:
                            replacement = self._expand(child, seen)
                    except (OSError, ValueError) as error:
                        line = source.count('\n', 0, directive.start()) + 1
                        raise ValueError(f'{path}:{line}: {error}') from error
                    end = directive.start('rest') + match.end()
            elif kind == 'pragma':
                match = ONCE.fullmatch(rest.rstrip())
                if match:
                    replacement = ''
                    end = directive.start('rest') + len(rest.rstrip())

            if replacement is not None:
                start = masked.index('#', directive.start(), directive.end())
                # 行頭・行末のコメントを残す。複数行コメントの境界も壊さない。
                output.append(source[cursor:start])
                output.append(replacement)
                cursor = end
        output.append(source[cursor:])
        body = ''.join(output)
        if root:
            return body
        # 条件付き include は別の分岐にも出現し得るため、コンパイラ側でも一度に制限。
        guard = f'CPPCPLIB_EXPANDER_FILE_{self.files[path]}_INCLUDED'
        return f'#ifndef {guard}\n#define {guard}\n{body}#endif // {guard}\n'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, help='展開する C++ ファイル')
    parser.add_argument('-o', '--output', type=Path, help='出力先（省略時は標準出力）')
    parser.add_argument('-I', '--include', action='append', default=[], metavar='DIR',
                        help='追加のヘッダ検索ディレクトリ（複数指定可）')
    args = parser.parse_args(argv)
    expander = Expander(args.include)
    try:
        result = expander.expand(args.source)
        if args.output is None or str(args.output) == '-':
            sys.stdout.write(result)
        else:
            if args.output.resolve() in expander.files:
                raise ValueError('出力先には入力ファイルや展開対象のヘッダを指定できません')
            args.output.write_text(result, encoding='utf-8')
    except (OSError, UnicodeError, ValueError) as error:
        print(f'expander.py: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
