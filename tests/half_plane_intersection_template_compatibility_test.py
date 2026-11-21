"""半平面交差と両テンプレートの併用を、include と using の順序を変えて確認する。"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import permutations
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
HEADERS = ["cplib/geometry/half_plane_convex_polygon.hpp"]
BODY = r"""
int main() {
    ll count = 3;
    vector<Int> values{1, 2, 3};
    rep(i, count) {
        assert(values[i] == i + 1);
    }

    vector<Line<double>> square{
        {{0, 0}, {2, 0}},
        {{2, 0}, {2, 2}},
        {{2, 2}, {0, 2}},
        {{0, 2}, {0, 0}}
    };
    HalfPlaneConvexPolygon<double> polygon(square);
    assert(polygon.contains({1, 1}, true));
    assert(polygon.on_edge({0, 1}));
    assert(polygon.half_planes().size() == 4);

    auto bounded = half_plane_intersection(square);
    assert(bounded.kind == hpiBounded);
    assert(bounded.dimension == 2);
    assert(bounded.vertices.size() == 4);
    assert(bounded.constraints.size() == square.size());
    assert(bounded.contains(Point<double>{1, 1}));
    assert(contains(bounded, Point<double>{0, 1}));
    assert(!bounded.contains(Point<double>{3, 1}));

    auto convex = bounded.toConvexPolygon();
    assert(convex.size() == 4);
    assert(convex.contains(Point<double>{1, 1}));
    assert(convex.diameter() == 8);
    assert(convex.toPolygon().size() == 4);
    assert(toConvexPolygon(bounded).size() == 4);

    HalfPlaneIntersection<double> empty;
    assert(empty.toConvexPolygon().size() == 0);
    assert(toConvexPolygon(empty).size() == 0);

    vector<Line<double>> half_plane{{{0, 0}, {1, 0}}};
    auto unbounded = half_plane_intersection(half_plane);
    assert(unbounded.kind == hpiUnbounded);
    assert(unbounded.dimension == 2);
    assert(unbounded.vertices.size() == 0);
    assert(unbounded.constraints.size() == 1);
    assert(contains(unbounded, Point<double>{100, 1}));
    assert(!unbounded.contains(Point<double>{0, -1}));

    bool member_threw = false;
    try {
        (void)unbounded.toConvexPolygon();
    } catch (const std::domain_error &) {
        member_threw = true;
    }
    assert(member_threw);

    bool free_function_threw = false;
    try {
        (void)toConvexPolygon(unbounded);
    } catch (const std::domain_error &) {
        free_function_threw = true;
    }
    assert(free_function_threw);

    using Rational = Fraction<Int>;
    Rational low(1, 3);
    Rational high(2, 3);
    vector<Line<Rational>> fractional_square{
        {{low, low}, {high, low}},
        {{high, low}, {high, high}},
        {{high, high}, {low, high}},
        {{low, high}, {low, low}}
    };
    auto exact = half_plane_intersection(fractional_square);
    assert(exact.kind == hpiBounded);
    assert(exact.dimension == 2);
    assert(exact.vertices.size() == 4);
    assert(exact.constraints.size() == fractional_square.size());
    assert(exact.contains(Point<Rational>{low, high}));
    assert(contains(exact, Point<Rational>{Rational(1, 2), Rational(1, 2)}));
    assert(!contains(exact, Point<Rational>{Rational(0), Rational(0)}));

    auto exact_convex = toConvexPolygon(exact);
    assert(exact_convex.size() == 4);
    assert(exact_convex.contains(Point<Rational>{Rational(1, 2), Rational(1, 2)}));
    assert(exact_convex.diameter() == Rational(2, 9));
    assert(exact_convex.toPolygon().size() == 4);
    assert(exact.toConvexPolygon().size() == 4);
}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="g++")
    parser.add_argument(
        "--syntax-only", action="store_true", help="実行せず構文検査だけを行う"
    )
    args = parser.parse_args()
    templates = ["cplib/tmpl/sheep.hpp", "cplib/tmpl/compact.hpp"]
    cases = [
        (order, library_first, using_first)
        for order in permutations(templates)
        for library_first in (False, True)
        for using_first in (False, True)
    ]

    with tempfile.TemporaryDirectory(prefix="half_plane_templates_") as directory:
        def check(index_case):
            index, (order, library_first, using_first) = index_case
            prefix = HEADERS if library_first else list(order)
            suffix = list(order) if library_first else HEADERS
            source = "".join(f"#include <{header}>\n" for header in prefix)
            if using_first:
                source += "using namespace cplib;\n"
            source += "".join(f"#include <{header}>\n" for header in suffix)
            if not using_first:
                source += "using namespace cplib;\n"
            source += BODY

            file = Path(directory) / f"case_{index}.cpp"
            binary = Path(directory) / f"case_{index}"
            file.write_text(source)
            command = [
                args.compiler, "-std=c++20", "-O1", "-I", str(ROOT / "src"), str(file)
            ]
            command += ["-fsyntax-only"] if args.syntax_only else ["-o", str(binary)]
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise AssertionError(
                    f"{order}, library_first={library_first}, "
                    f"using_first={using_first}\n{result.stderr}"
                )
            if not args.syntax_only:
                subprocess.run([str(binary)], check=True)

        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(check, enumerate(cases)))

    mode = "compile" if args.syntax_only else "runtime"
    print(f"{len(cases)} half-plane/template include-order {mode} checks passed")


if __name__ == "__main__":
    main()
