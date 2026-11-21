#include <cplib/geometry/geometry.hpp>
#include <random>

using namespace cplib;

using Integer = long long;
using PolygonType = HalfPlaneConvexPolygon<Integer>;
using Constraint = PolygonType::constraint_type;
using Rational = Fraction<Integer>;

void testPolygon() {
    const PolygonType polygon{{0, 0, 0}, {0, -1, 0}, {2, 2, -8}, {1, 1, -4},
                              {0, 1, -8}, {-1, 0, 0}, {0, -2, 0}, {1, -1, -4}};
    assert(polygon.dimension() == 2);
    assert(polygon.size() == 3);
    assert(polygon.contains({1, 1}, true));
    assert(polygon.contains({0, 1}));
    assert(!polygon.contains({0, 1}, true));
    assert(polygon.on_edge({0, 1}));
    assert(!polygon.contains({-1, 1}));
    assert(!polygon.on_edge({-1, 1}));
    assert(!polygon.contains({4, 1}));
    assert(area(toConvexPolygon(polygon).toPolygon()) == Rational(8));

    const auto left = convex_cut(polygon, Constraint{2, 0, -3});
    assert(left.size() == 4);
    assert(left.contains({1, 1}, true));
    assert(!left.contains({2, 1}));
    const Point<Rational> boundary{Rational(3, 2), Rational(1)};
    assert(contains(left, boundary));
    assert(!contains(left, boundary, true));
    assert(on_edge(left, boundary));
    assert(left.support_point({1, 0}).x == Rational(3, 2));
    assert(area(left.toPolygon()) == Rational(39, 8));

    const PolygonType right{{-2, 0, 3}, {1, 0, -4}, {0, -1, 0}, {0, 1, -4}};
    const auto segment = convex_intersection(left, right);
    assert(segment.dimension() == 1 && segment.size() == 2);
    assert(segment.half_planes().empty());
    assert(segment.contains(boundary) && !segment.contains(boundary, true));
    assert(!segment.contains({2, 1}));
    assert(segment.support_index({0, 1}) == 1);
    assert(segment.support_index({0, -1}) == 0);
    const auto point = convex_cut(segment, Constraint{0, 1, 0});
    assert(point.dimension() == 0 && point.size() == 1);
    assert(point.contains(point[0]) && !point.contains(point[0], true));
    assert(point.on_edge(point[0]));
    assert(convex_cut(point, Constraint{0, 1, 1}).empty());
    assert(convex_cut(polygon, Line<Integer>{{0, 0}, {0, 1}}).dimension() == 1);

    bool rejected = false;
    try {
        (void)PolygonType{{0, -1, 0}};
    } catch (const std::domain_error &) {
        rejected = true;
    }
    assert(rejected);
}

void testConversion() {
    for (const auto &vertices : {std::vector<Point<Integer>>{}, {{3, -2}}, {{-2, -3}, {1, 4}},
                                  {{1, -3}, {1, 4}}, {{-2, -3}, {1, -3}},
                                  {{0, 0}, {1, 0}, {1, 2}, {0, 2}}}) {
        const auto input = initConvexPolygon(vertices);
        const auto polygon = half_plane_convex_polygon(input);
        assert(polygon.size() == input.size());
        assert(std::distance(polygon.begin(), polygon.end()) == std::ssize(input));
        for (std::size_t i = 0; i < input.size(); ++i) {
            const Point<Rational> expected{input[i].x, input[i].y};
            assert(polygon[i] == expected);
            assert(polygon.contains(expected));
        }
    }

    const auto floating = half_plane_convex_polygon(initConvexPolygon(
        std::vector<Point<double>>{{-1, 0}, {2, 0}, {1, 2}, {-1, 2}}));
    assert(floating.size() == 4);
    assert(floating.contains({0, 1}, true));
    assert(floating.contains({0, 2}) && !floating.contains({0, 2}, true));
    assert(!floating.contains({2, 2}));
    const auto clipped = convex_cut(floating, Line<double>{{0, 1}, {1, 1}});
    assert(clipped.size() == 4);
    assert(clipped.contains({0, 1}) && clipped.on_edge({0, 1}));
    assert(!clipped.contains({0, 0}));
    assert(clipped.support_point({0, 1}).y == 2);

    const auto fractional = half_plane_convex_polygon(initConvexPolygon(
        std::vector<Point<Rational>>{{0, 0}, {Rational(2, 3), 0}, {0, Rational(3, 4)}}));
    assert(area(fractional.toConvexPolygon().toPolygon()) == Rational(1, 4));
    assert(fractional.contains({Rational(1, 4), Rational(1, 4)}, true));

    const HalfPlaneConvexPolygon<BigInt> big{{0, -1, 0}, {-1, 0, 0}, {1, 1, -2}};
    assert(big.size() == 3 && big.contains({1, 1}));
    assert(big.support_point({1, 0}).x == Fraction<BigInt>(2));
}

template <class T> void testIntegerType() {
    using P = HalfPlaneConvexPolygon<T>;
    const P polygon{{T(0), T(-1), T(0)}, {T(-1), T(0), T(0)}, {T(1), T(1), T(-2)}};
    assert(polygon.contains({T(0), T(0)}));
    assert(polygon.support_point({T(1), T(0)}).x == Fraction<T>(T(2)));

    const auto segment = half_plane_convex_polygon(ConvexPolygon<T>{{T(0), T(0)}, {T(1), T(0)}});
    assert(segment.dimension() == 1);
    assert(segment[1].x == Fraction<T>(T(1)));
}

void testRandomPolygons() {
    std::mt19937 random(49732);
    std::uniform_int_distribution<Integer> coefficient(-5, 5);
    for (int trial = 0; trial < 1500; ++trial) {
        std::vector<Constraint> constraints{{1, 0, -4}, {-1, 0, -4}, {0, 1, -4}, {0, -1, -4}};
        for (unsigned j = random() % 12; j > 0; --j)
            constraints.push_back({coefficient(random), coefficient(random), coefficient(random)});
        std::shuffle(constraints.begin(), constraints.end(), random);
        const auto expected = half_plane_intersection(constraints);
        const PolygonType polygon(constraints);
        assert(polygon.dimension() == expected.dimension);
        assert(polygon.size() == expected.vertices.size());
        assert((polygon.toConvexPolygon() == expected.toConvexPolygon()));

        if (polygon.dimension() == 2) {
            assert(polygon.half_planes().size() == polygon.size());
            auto evaluate = [](const Constraint &h, const Point<Rational> &p) {
                return Rational(h[0]) * p.x + Rational(h[1]) * p.y + Rational(h[2]);
            };
            for (std::size_t i = 0; i < polygon.size(); ++i) {
                const auto &h = polygon.half_planes()[i];
                assert(std::find(constraints.begin(), constraints.end(), h) != constraints.end());
                assert(evaluate(h, polygon[i]) == Rational(0));
                assert(evaluate(h, polygon[(i + 1) % polygon.size()]) == Rational(0));
                assert(evaluate(h, polygon[(i + 2) % polygon.size()]) < Rational(0));
            }
        }
        for (Integer x = -5; x <= 5; ++x) {
            for (Integer y = -5; y <= 5; ++y) {
                const Point<Integer> p{x, y};
                const bool inside = expected.contains(p);
                const bool strict = polygon.dimension() == 2 &&
                    std::ranges::all_of(constraints, [&](const Constraint &h) {
                        return (h[0] == 0 && h[1] == 0) || h[0] * x + h[1] * y + h[2] < 0;
                    });
                assert(polygon.contains(p) == inside);
                assert(polygon.contains(p, true) == strict);
                assert(polygon.on_edge(p) == (inside && !strict));
            }
        }
        if (polygon.empty())
            continue;
        for (int j = 0; j < 5; ++j) {
            const Point<Integer> direction{coefficient(random), coefficient(random)};
            const Point<Rational> d{direction.x, direction.y};
            const auto maximum = dot(d, polygon.support_point(direction));
            for (const auto &vertex : expected.vertices)
                assert(dot(d, vertex) <= maximum);
        }
    }
}

int main() {
    testPolygon();
    testConversion();
    testIntegerType<short>();
    testIntegerType<__int128>();
    testRandomPolygons();
}
