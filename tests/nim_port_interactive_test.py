import pathlib
import select
import subprocess
import tempfile

root = pathlib.Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='cppcplib_interactive_') as directory:
    source = pathlib.Path(directory) / 'interactive.cpp'
    executable = pathlib.Path(directory) / 'interactive'
    source.write_text('''#include <cplib/tmpl/sheep_interactive.hpp>
int main() {
    cplib::print("ready");
    auto first = (cplib::ii)();
    cplib::print(first + 1);
    auto word = (cplib::si)();
    cplib::print(word);
    auto values = cplib::lii(2);
    cplib::print(values[0] + values[1]);
    cplib::Yes();
}
''')
    subprocess.run(['g++', '-std=c++20', '-O2', '-I', str(root / 'src'),
                    str(source), '-o', str(executable)], check=True)
    process = subprocess.Popen([str(executable)], stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=0)

    def receive(expected):
        ready, _, _ = select.select([process.stdout], [], [], 5)
        assert ready, f'対話出力を待って停止しました: {expected!r}'
        actual = process.stdout.readline()
        assert actual == expected, (actual, expected)

    try:
        receive(b'ready\n')
        process.stdin.write(b'41\n')
        process.stdin.flush()
        receive(b'42\n')
        process.stdin.write(b'hello\n')
        process.stdin.flush()
        receive(b'hello\n')
        process.stdin.write(b'3 7\n')
        process.stdin.flush()
        receive(b'10\n')
        receive(b'Yes\n')
        assert process.wait(timeout=5) == 0, process.stderr.read()
    finally:
        if process.poll() is None:
            process.kill()
        process.communicate()
print('interactive pipe handshake passed')
