#include <cplib/geometry/half_plane_intersection.hpp>
#include <cplib/math/bigint.hpp>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string_view>

using namespace cplib;

namespace {

using Integer = long long;
using Rational = Fraction<__int128_t>;
using RationalPoint = Point<Rational>;

struct IntegerPoint {
    Integer x, y;
};

struct IntegerLine {
    IntegerPoint s, t;
};

struct OracleResult {
    HalfPlaneIntersectionKind kind = hpiEmpty;
    Int dimension = -1;
    std::vector<RationalPoint> vertices;
};

std::size_t checkedCases = 0;
const std::vector<IntegerLine> *currentInput = nullptr;

// 失敗した乱択入力を出力する。成功時 O(1)、失敗時 O(N)。
void require(bool condition, std::string_view message) {
    if (condition) {
        return;
    }
    std::cerr << "half-plane test case " << checkedCases << ": " << message << '\n';
    if (currentInput != nullptr) {
        for (const auto &line : *currentInput) {
            std::cerr << "((" << line.s.x << ", " << line.s.y << "), ("
                      << line.t.x << ", " << line.t.y << "))\n";
        }
    }
    std::abort();
}

// 整数点の差、外積、分数点の一致を厳密に計算する。各 O(1) 演算。
IntegerPoint direction(const IntegerLine &line) {
    return {line.t.x - line.s.x, line.t.y - line.s.y};
}

Integer integerCross(IntegerPoint a, IntegerPoint b) {
    return a.x * b.y - a.y * b.x;
}

bool samePoint(const RationalPoint &a, const RationalPoint &b) {
    return a.x == b.x && a.y == b.y;
}

// 与えられた整数制約をすべて満たすか直接判定する。O(N) 演算。
bool feasible(const std::vector<IntegerLine> &lines, const RationalPoint &point) {
    for (const auto &line : lines) {
        const auto d = direction(line);
        const Rational side = Rational(d.x) * (point.y - Rational(line.s.y))
                            - Rational(d.y) * (point.x - Rational(line.s.x));
        if (side < Rational(0)) {
            return false;
        }
    }
    return true;
}

// 全交点・原点・各境界への射影と recession 方向を列挙する独立解。O(N^3) 演算。
OracleResult oracle(const std::vector<IntegerLine> &lines) {
    OracleResult result;
    std::vector<RationalPoint> candidates{{Rational(0), Rational(0)}};
    std::vector<IntegerPoint> directions{{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto d = direction(lines[i]);
        directions.push_back(d);
        directions.push_back({-d.x, -d.y});
        const Integer c = integerCross(d, lines[i].s);
        const Integer squaredLength = d.x * d.x + d.y * d.y;
        candidates.push_back({Rational(-d.y * c, squaredLength),
                              Rational(d.x * c, squaredLength)});
        for (std::size_t j = 0; j < i; ++j) {
            const auto e = direction(lines[j]);
            const Integer determinant = integerCross(d, e);
            if (determinant == 0) {
                continue;
            }
            const Integer f = integerCross(e, lines[j].s);
            const RationalPoint point{Rational(c * e.x - d.x * f, determinant),
                                      Rational(e.y * c - d.y * f, determinant)};
            candidates.push_back(point);
            if (!feasible(lines, point)) {
                continue;
            }
            const bool duplicate = std::any_of(result.vertices.begin(), result.vertices.end(),
                [&](const RationalPoint &other) { return samePoint(point, other); });
            if (!duplicate) {
                result.vertices.push_back(point);
            }
        }
    }

    std::vector<RationalPoint> points;
    for (const auto &point : candidates) {
        if (feasible(lines, point)) {
            points.push_back(point);
        }
    }
    if (points.empty()) {
        return result;
    }

    std::vector<IntegerPoint> rays;
    for (const auto &candidate : directions) {
        bool valid = true;
        for (const auto &line : lines) {
            if (integerCross(direction(line), candidate) < 0) {
                valid = false;
                break;
            }
        }
        if (valid) {
            rays.push_back(candidate);
        }
    }
    result.kind = rays.empty() ? hpiBounded : hpiUnbounded;
    result.dimension = 0;
    std::vector<RationalPoint> spans;
    for (const auto &point : points) {
        spans.push_back({point.x - points.front().x, point.y - points.front().y});
    }
    for (const auto &ray : rays) {
        spans.push_back({Rational(ray.x), Rational(ray.y)});
    }
    for (const auto &a : spans) {
        if (a.x == Rational(0) && a.y == Rational(0)) {
            continue;
        }
        result.dimension = 1;
        for (const auto &b : spans) {
            if (a.x * b.y - a.y * b.x != Rational(0)) {
                result.dimension = 2;
                break;
            }
        }
        break;
    }
    return result;
}

// a*x+b*y>=c を整数端点の有向直線にする。O(1)。
IntegerLine constraint(Integer a, Integer b, Integer c) {
    if (b == 0) {
        require(a != 0 && c % a == 0, "invalid test constraint");
        const Integer x = c / a;
        if (a > 0) {
            return {{x, 1}, {x, 0}};
        }
        return {{x, 0}, {x, 1}};
    }
    require(b == 1 || b == -1, "unsupported test constraint");
    const Integer y = c / b;
    const Integer dy = -a / b;
    if (b > 0) {
        return {{0, y}, {1, y + dy}};
    }
    return {{1, y + dy}, {0, y}};
}

// 分数を浮動小数点へ変換して比較する。各 O(1)。
double asDouble(const Rational &value) {
    return static_cast<double>(value.num) / static_cast<double>(value.den);
}

bool near(double a, double b) {
    return std::abs(a - b) <= 1e-8 * std::max({1.0, std::abs(a), std::abs(b)});
}

// 並べ替えや重複で公開結果の頂点順序が変わらないことを調べる。O(N)。
template <class T>
void requireSameResult(const HalfPlaneIntersection<T> &a, const HalfPlaneIntersection<T> &b) {
    require(a.kind == b.kind && a.dimension == b.dimension, "order changed classification");
    require(a.vertices.size() == b.vertices.size(), "order changed vertex count");
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        require(a.vertices[i].x == b.vertices[i].x && a.vertices[i].y == b.vertices[i].y,
                "order changed canonical vertices");
    }
}

// 独立解との一致、包含判定、順序と重複への不変性を確認する。O(N^3) 演算。
void checkCase(const std::vector<IntegerLine> &lines, bool checkFloating = true) {
    ++checkedCases;
    currentInput = &lines;
    const auto expected = oracle(lines);
    std::vector<Line<Rational>> exactInput;
    std::vector<Line<double>> floatingInput;
    for (const auto &line : lines) {
        exactInput.push_back({{Rational(line.s.x), Rational(line.s.y)},
                              {Rational(line.t.x), Rational(line.t.y)}});
        floatingInput.push_back({{double(line.s.x), double(line.s.y)},
                                 {double(line.t.x), double(line.t.y)}});
    }
    const auto exact = half_plane_intersection(exactInput);
    require(exact.kind == expected.kind, "exact classification disagrees with oracle");
    require(exact.dimension == expected.dimension, "exact dimension disagrees with oracle");
    require(exact.vertices.size() == expected.vertices.size(), "exact vertex count differs");
    require(exact.constraints.size() == lines.size(), "exact constraints were discarded");
    ConvexPolygon<Rational> converted;
    if (exact.kind != hpiUnbounded) {
        converted = exact.toConvexPolygon();
        require(converted.size() == expected.vertices.size(), "conversion lost extreme points");
        if (!converted.empty()) {
            Rational expectedDiameter(0);
            for (const auto &a : expected.vertices) {
                for (const auto &b : expected.vertices) {
                    const Rational dx = a.x - b.x;
                    const Rational dy = a.y - b.y;
                    expectedDiameter = std::max(expectedDiameter, dx * dx + dy * dy);
                }
            }
            require(converted.diameter() == expectedDiameter, "converted diameter differs from oracle");
        }
    }
    for (const auto &point : exact.vertices) {
        const bool found = std::any_of(expected.vertices.begin(), expected.vertices.end(),
            [&](const RationalPoint &other) { return samePoint(point, other); });
        require(found, "unexpected exact vertex");
        require(exact.contains(point), "exact vertex is not contained");
    }
    for (Integer x = -3; x <= 3; ++x) {
        for (Integer y = -3; y <= 3; ++y) {
            const RationalPoint point{Rational(x), Rational(y)};
            const bool expectedContains = feasible(lines, point);
            require(exact.contains(point) == expectedContains, "incorrect exact contains");
            require(contains(exact, point) == expectedContains, "incorrect free exact contains");
            if (exact.kind != hpiUnbounded) {
                require(converted.contains(point) == expectedContains, "converted contains differs from oracle");
            }
        }
    }
    auto reordered = exactInput;
    std::reverse(reordered.begin(), reordered.end());
    requireSameResult(exact, half_plane_intersection(reordered));
    reordered.insert(reordered.end(), exactInput.begin(), exactInput.end());
    requireSameResult(exact, half_plane_intersection(reordered));
    std::mt19937 orderRng(static_cast<unsigned>(checkedCases));
    std::shuffle(reordered.begin(), reordered.end(), orderRng);
    requireSameResult(exact, half_plane_intersection(reordered));

    if (checkFloating) {
        const auto actual = half_plane_intersection(floatingInput);
        require(actual.kind == expected.kind, "floating classification disagrees with oracle");
        require(actual.dimension == expected.dimension, "floating dimension disagrees with oracle");
        require(actual.vertices.size() == expected.vertices.size(), "floating vertex count differs");
        require(actual.constraints.size() == lines.size(), "floating constraints were discarded");
        for (const auto &point : expected.vertices) {
            const bool found = std::any_of(actual.vertices.begin(), actual.vertices.end(),
                [&](const Point<double> &other) {
                    return near(asDouble(point.x), other.x) && near(asDouble(point.y), other.y);
                });
            require(found, "expected floating vertex is absent");
        }
        if (actual.kind == hpiBounded && actual.dimension == 2) {
            double doubledArea = 0;
            for (std::size_t i = 0; i < actual.vertices.size(); ++i) {
                const auto &a = actual.vertices[i];
                const auto &b = actual.vertices[(i + 1) % actual.vertices.size()];
                doubledArea += a.x * b.y - a.y * b.x;
            }
            require(doubledArea > 0, "bounded vertices are not counterclockwise");
        }
        auto floatingReordered = floatingInput;
        std::reverse(floatingReordered.begin(), floatingReordered.end());
        requireSameResult(actual, half_plane_intersection(floatingReordered));
        floatingReordered.insert(floatingReordered.end(), floatingInput.begin(), floatingInput.end());
        std::shuffle(floatingReordered.begin(), floatingReordered.end(), orderRng);
        requireSameResult(actual, half_plane_intersection(floatingReordered));
        for (Integer x = -3; x <= 3; ++x) {
            for (Integer y = -3; y <= 3; ++y) {
                const bool expectedContains = feasible(lines, {Rational(x), Rational(y)});
                const Point<double> point{double(x), double(y)};
                require(actual.contains(point) == expectedContains, "incorrect floating contains");
                require(contains(actual, point) == expectedContains, "incorrect free floating contains");
            }
        }
    }
    currentInput = nullptr;
}

// 空集合、全平面、半平面、帯、直線、半直線、線分、点、多角形を確認する。O(1) ケース。
void testBasicRegions() {
    checkCase({});
    checkCase({constraint(0, 1, 0)});
    checkCase({constraint(0, 1, 0), constraint(0, -1, -1)});
    checkCase({constraint(0, 1, 0), constraint(0, -1, 0)});
    checkCase({constraint(0, 1, 1), constraint(0, -1, 0)});
    checkCase({constraint(1, 0, 0), constraint(-1, 0, 0)});
    checkCase({constraint(1, 0, 0), constraint(-1, 0, 0), constraint(0, 1, 0)});
    checkCase({constraint(1, 0, 0), constraint(-1, 0, 0),
               constraint(0, 1, 0), constraint(0, -1, -2)});
    checkCase({constraint(1, 0, 0), constraint(-1, 0, 0),
               constraint(0, 1, 0), constraint(0, -1, 0)});
    checkCase({constraint(0, 1, 0), constraint(0, -1, 0), constraint(1, 0, 0)});
    checkCase({constraint(0, 1, 0), constraint(0, -1, 0),
               constraint(1, 0, 0), constraint(-1, 0, -2)});
    checkCase({constraint(-1, 1, 0), constraint(1, -1, 0),
               constraint(1, 0, -2), constraint(-1, 0, -3)});
    checkCase({constraint(1, 1, 0), constraint(-1, 1, 0), constraint(0, -1, 0)});
    checkCase({constraint(1, 1, 1), constraint(-1, 1, 1), constraint(0, -1, 0)});
    checkCase({constraint(1, 0, 0), constraint(-1, 0, -3),
               constraint(0, 1, 0), constraint(0, -1, -2)});
    checkCase({constraint(0, 1, 0), constraint(0, 1, 0),
               constraint(0, 1, 1), constraint(0, 1, -1)});

    // 非二進分数の境界でできる点は、浮動小数点の退化分類を要求しない。
    checkCase({{{0, 0}, {-1, 2}}, {{-2, -3}, {0, 2}}, {{-1, -2}, {-3, 3}},
               {{-3, 1}, {2, 1}}, {{-2, -2}, {2, -3}}, {{2, 0}, {-3, 3}},
               {{2, -2}, {3, -3}}, {{-2, 0}, {-2, 3}}}, false);
}

// 任意端点と構造的な整数係数を使って独立解と比較する。O(K N^3) 演算。
void testRandomRegions() {
    std::mt19937_64 rng(20261002);
    auto randomInteger = [&](Integer minimum, Integer maximum) {
        return minimum + static_cast<Integer>(rng() % (maximum - minimum + 1));
    };
    for (Int iteration = 0; iteration < 2000; ++iteration) {
        std::vector<IntegerLine> lines;
        const Integer count = randomInteger(0, 10);
        for (Integer i = 0; i < count; ++i) {
            const Integer a = randomInteger(-3, 3);
            const Integer b = randomInteger(-1, 1);
            if (a == 0 && b == 0) {
                continue;
            }
            const Integer c = randomInteger(-4, 4) * (b == 0 ? std::abs(a) : 1);
            lines.push_back(constraint(a, b, c));
            if (randomInteger(0, 4) == 0) {
                lines.push_back(lines.back());
            }
        }
        checkCase(lines);
    }
    for (Int iteration = 0; iteration < 2000; ++iteration) {
        std::vector<IntegerLine> lines;
        const Integer count = randomInteger(0, 10);
        for (Integer i = 0; i < count; ++i) {
            const IntegerPoint s{randomInteger(-3, 3), randomInteger(-3, 3)};
            const IntegerPoint t{randomInteger(-3, 3), randomInteger(-3, 3)};
            if (s.x != t.x || s.y != t.y) {
                lines.push_back({s, t});
            }
        }
        checkCase(lines, false);
    }
    for (Int iteration = 0; iteration < 400; ++iteration) {
        std::vector<IntegerLine> lines{constraint(1, 0, -8), constraint(-1, 0, -8),
                                       constraint(0, 1, -8), constraint(0, -1, -8)};
        const Integer count = randomInteger(0, 8);
        for (Integer i = 0; i < count; ++i) {
            const Integer a = randomInteger(-3, 3);
            const Integer b = randomInteger(-1, 1);
            if (a == 0 && b == 0) {
                continue;
            }
            const Integer c = randomInteger(-4, 4) * (b == 0 ? std::abs(a) : 1);
            lines.push_back(constraint(a, b, c));
        }
        checkCase(lines);
    }
}

// 各浮動小数点型と通常の整数を基底とする分数でも同じAPIを使えることを確認する。O(1)。
template <class T>
void testCoordinateType() {
    const std::vector<Line<T>> square{
        {{T(0), T(0)}, {T(1), T(0)}}, {{T(1), T(0)}, {T(1), T(1)}},
        {{T(1), T(1)}, {T(0), T(1)}}, {{T(0), T(1)}, {T(0), T(0)}}};
    const auto result = half_plane_intersection(square);
    require(result.kind == hpiBounded && result.dimension == 2, "coordinate type failed");
    require(result.vertices.size() == 4, "coordinate type lost square vertices");
    require(result.contains(Point<T>{T(0), T(0)}), "closed boundary is excluded");
    require(!contains(result, Point<T>{T(2), T(0)}), "outside square was included");
}

// 変換後の凸多角形の操作、退化領域、コピーの独立性と非有界領域の拒否を調べる。O(1)。
template <class T>
void testConvexPolygonConversion() {
    const Line<T> bottom{{T(0), T(0)}, {T(1), T(0)}};
    const Line<T> top{{T(1), T(2)}, {T(0), T(2)}};
    const Line<T> left{{T(0), T(1)}, {T(0), T(0)}};
    const Line<T> right{{T(2), T(0)}, {T(2), T(1)}};
    const Line<T> reverseBottom{bottom.t, bottom.s};
    const Line<T> reverseLeft{left.t, left.s};
    const auto square = half_plane_intersection<T>({bottom, right, top, left});
    auto polygon = square.toConvexPolygon();
    static_assert(std::is_same_v<decltype(polygon), ConvexPolygon<T>>);
    const auto other = toConvexPolygon(square);
    require(polygon.size() == 4 && other.size() == 4, "square conversion lost vertices");
    require(polygon.contains(Point<T>{T(1), T(1)}, true), "converted square lost its interior");
    require(polygon.on_edge(Point<T>{T(0), T(1)}), "converted square lost its boundary");
    require(!polygon.contains(Point<T>{T(3), T(1)}), "converted square included an outside point");
    require(polygon.diameter() == T(8), "converted square has wrong squared diameter");
    require(polygon.toPolygon().size() == 4, "conversion cannot be chained to Polygon");
    polygon[0] = {T(9), T(9)};
    require(exact_equal(square.vertices[0], Point<T>{T(0), T(0)}), "conversion changed the source");
    require(exact_equal(other[0], Point<T>{T(0), T(0)}), "converted polygons share storage");

    const auto empty = half_plane_intersection<T>({{top.t, top.s}, reverseBottom});
    require(empty.kind == hpiEmpty && toConvexPolygon(empty).empty(), "empty region conversion failed");
    require(HalfPlaneIntersection<T>{}.toConvexPolygon().empty(), "default empty conversion failed");

    const auto point = half_plane_intersection<T>({bottom, reverseBottom, left, reverseLeft});
    const auto convertedPoint = point.toConvexPolygon();
    require(convertedPoint.size() == 1 && convertedPoint.diameter() == T(0), "point conversion failed");
    require(convertedPoint.contains(Point<T>{T(0), T(0)}), "converted point is excluded");
    require(!convertedPoint.contains(Point<T>{T(0), T(1)}), "converted point includes another point");

    const Line<T> diagonal{{T(0), T(0)}, {T(1), T(1)}};
    const auto segment = half_plane_intersection<T>({diagonal, {diagonal.t, diagonal.s}, left, right});
    const auto convertedSegment = toConvexPolygon(segment);
    require(convertedSegment.size() == 2 && convertedSegment.diameter() == T(8), "segment conversion failed");
    require(convertedSegment.contains(Point<T>{T(1), T(1)}), "converted segment lost its interior");
    require(!convertedSegment.contains(Point<T>{T(1), T(0)}), "converted segment includes an off-line point");
    const auto verticalSegment = half_plane_intersection<T>({left, reverseLeft, bottom, top});
    const auto convertedVertical = toConvexPolygon(verticalSegment);
    require(convertedVertical.size() == 2 && convertedVertical.diameter() == T(4), "vertical segment conversion failed");

    const std::vector<std::vector<Line<T>>> unboundedInputs{
        {}, {bottom}, {bottom, top}, {bottom, reverseBottom},
        {bottom, reverseBottom, left}, {bottom, left}};
    for (const auto &lines : unboundedInputs) {
        const auto region = half_plane_intersection(lines);
        require(region.kind == hpiUnbounded, "unbounded conversion test has bounded input");
        bool memberRejected = false;
        bool freeRejected = false;
        try {
            (void)region.toConvexPolygon();
        } catch (const std::domain_error &) {
            memberRejected = true;
        }
        try {
            (void)toConvexPolygon(region);
        } catch (const std::domain_error &) {
            freeRejected = true;
        }
        require(memberRejected && freeRejected, "unbounded region was silently converted");
    }
}

// 面積の演算がアンダーフローする大きさでも、変換は頂点を保つ。O(1)。
void testTinyPolygonConversion() {
    const double tiny = 1e-200;
    const auto region = half_plane_intersection<double>({
        {{0, 0}, {1, 0}}, {{tiny, 0}, {tiny, 1}},
        {{1, tiny}, {0, tiny}}, {{0, 1}, {0, 0}}});
    require(region.kind == hpiBounded && region.dimension == 2, "tiny square construction failed");
    const auto polygon = region.toConvexPolygon();
    require(polygon.size() == 4, "conversion collapsed a tiny square");
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        require(exact_equal(polygon[i], region.vertices[i]), "conversion changed a tiny coordinate");
    }
}

// 配列、初期化リスト、const span の入力を確認する。O(1) サイズの制約。
void testInputContainers() {
    const std::array<Line<double>, 4> square{{
        {{0, 0}, {1, 0}}, {{1, 0}, {1, 1}},
        {{1, 1}, {0, 1}}, {{0, 1}, {0, 0}}}};
    const auto fromArray = half_plane_intersection(square);
    const std::span<const Line<double>> view(square.data(), square.size());
    const auto fromSpan = half_plane_intersection(view);
    const auto fromList = half_plane_intersection<double>({
        {{0, 0}, {1, 0}}, {{1, 0}, {1, 1}},
        {{1, 1}, {0, 1}}, {{0, 1}, {0, 0}}});
    require(fromArray.kind == hpiBounded && fromArray.dimension == 2,
            "array input failed");
    require(fromArray.vertices.size() == 4, "array input lost vertices");
    requireSameResult(fromArray, fromSpan);
    requireSameResult(fromArray, fromList);
}

// 真の分数座標と、EPSより狭い領域・ほぼ平行な境界を確認する。O(1) ケース。
void testFractionalAndNarrowRegions() {
    const Rational oneThird(1, 3);
    const Rational twoThirds(2, 3);
    const std::vector<RationalPoint> corners{{oneThird, oneThird}, {twoThirds, oneThird},
                                             {twoThirds, twoThirds}, {oneThird, twoThirds}};
    std::vector<Line<Rational>> lines;
    for (std::size_t i = 0; i < corners.size(); ++i) {
        lines.push_back({corners[i], corners[(i + 1) % corners.size()]});
    }
    const auto square = half_plane_intersection(lines);
    require(square.kind == hpiBounded && square.dimension == 2, "fractional square failed");
    require(square.vertices.size() == 4, "fractional square lost vertices");
    for (const auto &point : square.vertices) {
        require((point.x == oneThird || point.x == twoThirds)
                && (point.y == oneThird || point.y == twoThirds), "wrong fractional vertex");
        require(square.contains(point), "fractional boundary is excluded");
    }

    const double tiny = 1e-12;
    const std::vector<Line<double>> narrowInput{
        {{0, 0}, {1, 0}}, {{1, tiny}, {0, tiny}}};
    const std::vector<Line<double>> nearParallelInput{
        {{0, 0}, {1, 0}}, {{1, 1 - tiny}, {0, 1}}};
    const auto referenceNarrow = half_plane_intersection(narrowInput);
    const auto referenceParallel = half_plane_intersection(nearParallelInput);
    const long double savedEps = GEOMETRY_EPS;
    for (const long double eps : {0.0L, 1e-10L, 1.0L, 1e12L}) {
        GEOMETRY_EPS = eps;
        const auto narrow = half_plane_intersection(narrowInput);
        require(narrow.kind == hpiUnbounded && narrow.dimension == 2, "EPS collapsed a narrow strip");
        require(narrow.vertices.empty(), "strip has unexpected extreme points");
        require(narrow.contains(Point<double>{0, tiny / 2}), "narrow strip lost its interior");
        require(!narrow.contains(Point<double>{0, -tiny}), "EPS admitted an outside point");
        require(!narrow.contains(Point<double>{0, 2 * tiny}), "EPS enlarged the strip");
        requireSameResult(referenceNarrow, narrow);
        const auto nearParallel = half_plane_intersection(nearParallelInput);
        require(nearParallel.kind == hpiUnbounded && nearParallel.dimension == 2,
                "EPS treated distinct slopes as parallel");
        require(nearParallel.vertices.size() == 1 && nearParallel.vertices.front().x > 1e11,
                "near-parallel intersection is absent");
        requireSameResult(referenceParallel, nearParallel);
    }
    GEOMETRY_EPS = savedEps;
}

// 無効な座標、同一点の端点、非有限の中間値を例外として拒否する。O(1) ケース。
void testInvalidInput() {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<Line<double>> invalidLines{
        {{0, 0}, {0, 0}}, {{nan, 0}, {1, 1}}, {{0, 0}, {infinity, 1}},
        {{-infinity, 0}, {0, 1}}, {{-1e308, 0}, {1e308, 1}},
        {{0, 0}, {1e-308, 1e308}}, {{1e200, 1e200}, {2e200, 0}}};
    for (const auto &line : invalidLines) {
        bool rejected = false;
        try {
            (void)half_plane_intersection(std::vector<Line<double>>{line});
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "invalid floating line was accepted");
    }
    const std::vector<Rational> invalidFractions{
        Rational(), Rational(1, 0, false), Rational(1, -1, false)};
    for (const auto &bad : invalidFractions) {
        bool rejected = false;
        try {
            const Line<Rational> line{{bad, Rational(0)}, {Rational(1), Rational(1)}};
            (void)half_plane_intersection(std::vector<Line<Rational>>{line});
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "invalid fraction was accepted");
    }
    const std::vector<HalfPlaneIntersection<double>> regions{
        half_plane_intersection(std::vector<Line<double>>{}),
        half_plane_intersection(std::vector<Line<double>>{
            {{0, 1}, {1, 1}}, {{1, 0}, {0, 0}}})};
    for (const auto &region : regions) {
        for (const Point<double> point : {Point<double>{nan, 0}, Point<double>{0, infinity}}) {
            bool rejected = false;
            try {
                (void)region.contains(point);
            } catch (const std::invalid_argument &) {
                rejected = true;
            }
            require(rejected, "contains accepted a nonfinite point");
        }
    }
    const auto plane = half_plane_intersection(std::vector<Line<Rational>>{});
    for (const auto &bad : invalidFractions) {
        bool rejected = false;
        try {
            (void)contains(plane, RationalPoint{bad, Rational(0)});
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "contains accepted an invalid fraction");
    }
}

// 同方向の大量制約と大量の有効な極点を持つ入力を確認する。O(N log N) を想定。
void testLargeInputs() {
    constexpr Int count = 20000;
    std::vector<Line<double>> parallel;
    std::vector<Line<double>> staircase;
    parallel.reserve(count);
    staircase.reserve(count);
    for (Int i = 0; i < count; ++i) {
        const double k = static_cast<double>(i);
        parallel.push_back({{0, k}, {1, k}});
        staircase.push_back({{0, -k * k}, {1, k - k * k}});
    }
    const auto strip = half_plane_intersection(parallel);
    require(strip.kind == hpiUnbounded && strip.dimension == 2, "large parallel input failed");
    require(strip.vertices.empty() && strip.constraints.size() == count,
            "large parallel input lost constraints or gained vertices");
    require(strip.contains(Point<double>{0, count - 1.0}), "strongest parallel bound was lost");
    require(!strip.contains(Point<double>{0, count - 2.0}), "weaker parallel bound was used");

    const auto stairs = half_plane_intersection(staircase);
    require(stairs.kind == hpiUnbounded && stairs.dimension == 2, "large staircase failed");
    require(stairs.vertices.size() == count - 1, "large staircase lost active constraints");
    require(stairs.constraints.size() == count, "large staircase discarded input constraints");
    for (std::size_t i = 0; i < stairs.vertices.size(); ++i) {
        const double k = static_cast<double>(i);
        require(stairs.vertices[i].x == 2 * k + 1 && stairs.vertices[i].y == k * (k + 1),
                "large staircase has an incorrect vertex or vertex order");
    }
    for (const std::size_t index : {std::size_t(0), stairs.vertices.size() / 2,
                                   stairs.vertices.size() - 1}) {
        require(stairs.contains(stairs.vertices[index]), "large staircase vertex is infeasible");
    }
    std::reverse(staircase.begin(), staircase.end());
    requireSameResult(stairs, half_plane_intersection(staircase));
}

} // namespace

int main() {
    testBasicRegions();
    testRandomRegions();
    testCoordinateType<float>();
    testCoordinateType<double>();
    testCoordinateType<long double>();
    testCoordinateType<Fraction<Int>>();
    testCoordinateType<Rational>();
    testCoordinateType<Fraction<BigInt>>();
    testConvexPolygonConversion<float>();
    testConvexPolygonConversion<double>();
    testConvexPolygonConversion<long double>();
    testConvexPolygonConversion<Fraction<Int>>();
    testConvexPolygonConversion<Rational>();
    testConvexPolygonConversion<Fraction<BigInt>>();
    testTinyPolygonConversion();
    testInputContainers();
    testFractionalAndNarrowRegions();
    testInvalidInput();
    testLargeInputs();
    std::cout << "half-plane intersection: " << checkedCases << " oracle cases passed\n";
}
