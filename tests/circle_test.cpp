#include <cplib/geometry/circle.hpp>
#include <cplib/geometry/minimum_enclosing_circle.hpp>
#include <cplib/math/bigint.hpp>
#include <cplib/math/int128.hpp>
#include <random>
#include <optional>
#include <type_traits>
using namespace cplib;

bool close(double a, double b, double scale = 1.0) {
    return std::abs(a - b) <= 2e-8 * scale;
}
void checkTangents(const Circle<double> &a, const Circle<double> &b, Int expected, double scale = 1.0) {
    auto result = common_tangents(a, b);
    assert(result.kind == circleFinite && Int(result.tangents.size()) == expected);
    assert(common_tangent_count(a, b) == expected);
    auto reverse = common_tangents(b, a);
    assert(reverse.tangents.size() == result.tangents.size());
    for (const auto &line : result.tangents) {
        assert(close(std::hypot(line.first.x - a.center().x, line.first.y - a.center().y), a.radius(), scale));
        assert(close(std::hypot(line.second.x - b.center().x, line.second.y - b.center().y), b.radius(), scale));
        assert(close(std::hypot(line.direction.x, line.direction.y), 1));
        assert(std::abs(dot(line.first - a.center(), line.direction)) <= 2e-8 * scale);
        assert(std::abs(dot(line.second - b.center(), line.direction)) <= 2e-8 * scale);
        assert(std::abs(cross(line.second - line.first, line.direction)) <= 2e-8 * scale);
        bool found = false;
        for (const auto &other : reverse.tangents)
            if (std::hypot(line.first.x - other.second.x, line.first.y - other.second.y) <= 2e-8 * scale &&
                std::hypot(line.second.x - other.first.x, line.second.y - other.first.y) <= 2e-8 * scale)
                found = true;
        assert(found);
    }
}
void testFloating() {
    auto unit = initCircle(Point<double>{0, 0}, 1);
    auto zero = initCircle(Point<double>{0, 0}, 0);
    assert(unit.contains(Point<double>{1, 0}) && !unit.contains(Point<double>{2, 0}));
    assert(cross_points(unit, unit).kind == circleInfinite);
    assert(common_tangents(unit, unit).kind == circleInfinite);
    assert(tangent_count(zero, Point<int>{0, 0}) == -1);
    assert(cross_points(zero, zero).points.size() == 1);
    assert(cross_points(unit, initCircle(Point<double>{0, 0}, 1 + 1e-12)).points.empty());
    checkTangents(zero, initCircle(Point<double>{2, 0}, 0), 1);
    checkTangents(unit, initCircle(Point<double>{3, 0}, 1), 4);
    checkTangents(unit, initCircle(Point<double>{2, 0}, 1), 3);
    checkTangents(unit, initCircle(Point<double>{1, 0}, 1), 2);
    checkTangents(initCircle(Point<double>{0, 0}, 2), initCircle(Point<double>{1, 0}, 1), 1);
    checkTangents(initCircle(Point<double>{0, 0}, 2), initCircle(Point<double>{0.5, 0}, 1), 0);
    for (double scale : {1e-150, 1.0, 1e150}) {
        auto a = initCircle(Point<double>{0, 0}, scale);
        auto b = initCircle(Point<double>{1.2 * scale, 0}, scale);
        auto result = cross_points(a, b);
        assert(result.points.size() == 2);
        for (auto p : result.points) {
            assert(close(std::hypot(p.x, p.y), scale, scale));
            assert(close(std::hypot(p.x - 1.2 * scale, p.y), scale, scale));
        }
        checkTangents(a, b, 2, scale);
        assert(cross_points(a, Line<double>{{-2 * scale, 0}, {2 * scale, 0}}).points.size() == 2);
        assert(cross_points(a, Segment<double>{{-scale, 0}, {0, 0}}).points.size() == 1);
        assert(cross_points(a, Segment<double>{{scale, 0}, {scale, 0}}).points.size() == 1);
    }
    assert(cross_points(unit, initCircle(Point<double>{1e-200, 0}, 1)).points.size() == 2);
    for (double radius : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        bool rejected = false;
        try { initCircle(Point<double>{0, 0}, radius); }
        catch (const std::invalid_argument &) { rejected = true; }
        assert(rejected);
    }
}

template <class T> void testExactType() {
    Point<T> o{T(0), T(0)}, p{T(1), T(0)};
    auto circle = initCircle(o, p), diameter = initDiameterCircle(o, p);
    static_assert(std::is_same_v<decltype(circle), Circle<T>>);
    using B = detail::CircleBase<T>;
    static_assert(std::is_same_v<decltype(circle.center_exact()), Point<Fraction<B>>>);
    assert(circle.contains(o) && circle.on_circle(p));
    assert(circle.radius() == 1 && diameter.radius() == 0.5);
    assert(diameter.center_approx().x == 0.5);
    assert(intersection_count(circle, diameter) == 1);
    assert(common_tangent_count(circle, diameter) == 1);
    assert(tangent_count(diameter, p) == 1);
    assert(intersection_count(diameter, Line<T>{o, p}) == 2);
    assert(intersection_count(diameter, Segment<T>{o, p}) == 2);
    assert(cross_points_approx(diameter, circle).points.size() == 1);
    assert(tangent_lines_approx(circle, p).tangents.size() == 1);
    for (auto point : diameter.points_exact())
        assert(diameter.on_circle(point));
    assert(circle == initCircle(Point<T>{T(1), T(0)}, Point<T>{T(0), T(1)}, Point<T>{T(-1), T(0)}));
    assert(circle != diameter);
    auto enclosing = minimum_enclosing_circle(std::vector<Point<T>>{o, {T(4), T(0)}, {T(2), T(4)}});
    auto squared = enclosing.radius_squared_exact();
    assert(squared.num == B(25) && squared.den == B(4));
    for (auto point : enclosing.points_exact())
        assert(enclosing.on_circle(point));
    bool rejected = false;
    try { initCircle(o, p, Point<T>{T(2), T(0)}); }
    catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
}
void testRandomCircles(std::mt19937_64 &rng) {
    for (Int trial = 0; trial < 2000; ++trial) {
        Int x = Int(rng() % 41) - 20, y = Int(rng() % 41) - 20;
        Int ra = rng() % 11, rb = rng() % 11;
        auto a = initCircle(Point<Int>{0, 0}, ra), b = initCircle(Point<Int>{x, y}, rb);
        auto af = toFloatCircle(a), bf = toFloatCircle(b);
        Int distance = x * x + y * y, sum = (ra + rb) * (ra + rb), difference = (ra - rb) * (ra - rb);
        Int expected = 0;
        if (distance == 0) {
            if (ra == rb)
                expected = ra == 0 ? 1 : -1;
        } else if (distance == sum || distance == difference) {
            expected = 1;
        } else if (difference < distance && distance < sum) {
            expected = 2;
        }
        assert(intersection_count(a, b) == expected);
        auto result = cross_points(a, b);
        assert(expected == -1 ? result.kind == circleInfinite : Int(result.points.size()) == expected);
        for (auto p : result.points) {
            assert(close(std::hypot(p.x, p.y), double(ra)));
            assert(close(std::hypot(p.x - x, p.y - y), double(rb)));
        }
        Int tangents = common_tangent_count(a, b);
        auto lines = common_tangents(a, b);
        assert(tangents == -1 ? lines.kind == circleInfinite : Int(lines.tangents.size()) == tangents);
        assert(tangents == common_tangent_count(af, bf));
        if (tangents >= 0)
            checkTangents(af, bf, tangents);
        Line<Int> line{{-30, y}, {30, y}};
        Int count = y * y == ra * ra ? 1 : y * y < ra * ra ? 2 : 0;
        assert(intersection_count(a, line) == count);
        assert(Int(cross_points(a, line).points.size()) == count);
        Segment<Int> segment{{x, y}, {x + 1, y}};
        Int segmentCount = 0;
        for (auto p : cross_points(af, Line<double>{{-30, double(y)}, {30, double(y)}}).points)
            segmentCount += p.x >= x - 1e-9 && p.x <= x + 1 + 1e-9;
        assert(intersection_count(a, segment) == segmentCount);
        assert(Int(cross_points(a, segment).points.size()) == segmentCount);
        Point<Int> point{x, y};
        auto tangentsFromPoint = tangent_lines(a, point);
        Int pointCount = tangent_count(a, point);
        assert(pointCount == -1 ? tangentsFromPoint.kind == circleInfinite
                               : Int(tangentsFromPoint.tangents.size()) == pointCount);
    }
}
void testMinimum(std::mt19937_64 &rng) {
    using P = Point<Int128>;
    for (Int trial = 0; trial < 300; ++trial) {
        std::vector<P> points(1 + rng() % 9);
        for (auto &p : points)
            p = {Int128(Int(rng() % 9) - 4), Int128(Int(rng() % 9) - 4)};
        std::optional<Circle<Int128>> expected;
        auto consider = [&](const Circle<Int128> &candidate) {
            for (auto p : points)
                if (!candidate.contains(p))
                    return;
            if (!expected || detail::mecRadiusLess(candidate, *expected))
                expected = candidate;
        };
        for (std::size_t i = 0; i < points.size(); ++i) {
            consider(initCircle(points[i], points[i], points[i]));
            for (std::size_t j = 0; j < i; ++j) {
                consider(initDiameterCircle(points[i], points[j]));
                for (std::size_t k = 0; k < j; ++k) {
                    try { consider(initCircle(points[i], points[j], points[k])); }
                    catch (const std::invalid_argument &) {}
                }
            }
        }
        assert(expected);
        std::vector<Point<double>> floating;
        for (auto p : points)
            floating.push_back({double(p.x), double(p.y)});
        for (Int seed : {0, 1, -1, 111}) {
            auto actual = minimum_enclosing_circle(points, seed);
            assert(actual == *expected);
            auto approx = minimum_enclosing_circle(floating, seed);
            assert(close(approx.radius(), expected->radius_approx()));
            for (auto p : floating)
                assert(approx.contains(p));
        }
    }
    BigInt giant("1" + std::string(300, '0'));
    auto exact = minimum_enclosing_circle(std::vector<Point<BigInt>>{{giant, giant}, {giant + 1, giant}});
    assert(exact.radius_squared_exact() == Fraction<BigInt>(BigInt(1), BigInt(4)));
    using F = Fraction<BigInt>;
    F tiny(BigInt(1), giant);
    auto tinyCircle = minimum_enclosing_circle(std::vector<Point<F>>{{F(0), F(0)}, {tiny, F(0)}});
    auto squared = tinyCircle.radius_squared_exact();
    assert(squared.num == BigInt(1) && squared.den == giant * giant * 4);
    auto mixed = initCircle(Point<Fraction<Int>>{Fraction<Int>(-1, -2, false), Fraction<Int>(1, -2, false)}, 0.0);
    assert(mixed.center().x == 0.5 && mixed.center().y == -0.5);
}
int main() {
    std::mt19937_64 rng(20261115);
    testFloating();
    testExactType<Int>();
    testExactType<Int128>();
    testExactType<BigInt>();
    testExactType<Fraction<Int>>();
    testExactType<Fraction<Int128>>();
    testExactType<Fraction<BigInt>>();
    testRandomCircles(rng);
    testMinimum(rng);
}
