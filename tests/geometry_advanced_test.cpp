#include <cplib/geometry/geometry.hpp>
#include <cmath>
#include <iostream>
#include <random>
#include <set>

using namespace cplib;

namespace {
using P = Point<Int>;
using R = Fraction<BigInt>;
using RP = Point<R>;

Int random_coordinate(std::mt19937 &random, int limit = 20) {
    return Int(random() % (2 * limit + 1)) - limit;
}

std::vector<P> random_points(std::mt19937 &random, std::size_t n, int limit = 20) {
    std::vector<P> points(n);
    for (auto &p : points)
        p = {random_coordinate(random, limit), random_coordinate(random, limit)};
    return points;
}

auto hull(const std::vector<P> &points) { return initConvexPolygon(convex_hull(points)); }

void test_convex_operations(std::mt19937 &random) {
    const auto square = hull({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
    const auto triangle = hull({{0, 0}, {4, 0}, {4, 4}});
    const auto diagonal = hull({{-1, 5}, {5, -1}});
    const auto clipped = convex_intersection(square, diagonal);
    assert((clipped.size() == 2 && clipped.front() == RP{R(0), R(4)} && clipped.back() == RP{R(4), R(0)}));
    assert(convex_intersection_area(square, triangle) == R(8));
    assert(convex_intersection(square, hull({{4, 4}, {5, 4}, {5, 5}})).size() == 1);
    assert(convex_intersection(square, hull({{5, 5}})).empty());
    assert(convex_intersection(square, hull({{3, 3}})).size() == 1);
    assert(convex_intersection(hull({{-1, 0}, {1, 0}}), hull({{0, -1}, {0, 1}})).size() == 1);
    assert(convex_intersection_area(hull({{0, 0}, {2, 0}, {0, 2}}), hull({{1, 0}, {3, 0}, {1, 2}})) == R(1, 2));
    assert(convex_intersection_area(ConvexPolygon<double>{{0, 0}, {2, 0}, {0, 2}},
                                    ConvexPolygon<double>{{1, 0}, {3, 0}, {1, 2}}) == 0.5);

    const std::vector<P> seeds{{0, 0}, {4, 0}, {4, 4}, {0, 4}, {2, 0}, {4, 2}};
    for (auto p : seeds)
        assert(!convex_tangents(square, p));
    auto shared_edge = convex_line_intersection(square, Line<Int>{{5, 4}, {-1, 4}});
    assert((shared_edge && shared_edge->s == RP{R(4), R(4)} && shared_edge->t == RP{R(0), R(4)}));
    assert(!convex_line_intersection(square, Line<Int>{{5, 5}, {-1, 5}}));

    for (int round = 0; round < 800; ++round) {
        auto a = hull(random_points(random, random() % 30));
        auto b = hull(random_points(random, random() % 30));
        std::vector<P> sums;
        for (auto p : a)
            for (auto q : b)
                sums.push_back(p + q);
        assert(minkowski_sum(a, b) == hull(sums));

        for (auto p : random_points(random, 20, 24)) {
            auto tangents = convex_tangents(a, p);
            assert(tangents.has_value() == (!a.empty() && !a.contains(p)));
            if (tangents)
                for (auto q : a) {
                    assert(cross(a[tangents->first] - p, q - p) >= 0);
                    assert(cross(a[tangents->second] - p, q - p) <= 0);
                }

            P direction{random_coordinate(random, 4), random_coordinate(random, 4)};
            if (direction == P{0, 0})
                direction.x = 1;
            auto intersection = convex_line_intersection(a, Line<Int>{p, p + direction});
            if (a.empty()) {
                assert(!intersection);
                continue;
            }
            auto rational = [](P point) -> RP { return {R(point.x), R(point.y)}; };
            std::vector<RP> expected;
            for (std::size_t i = 0; i < a.size(); ++i) {
                P u = a[i], v = a[(i + 1) % a.size()];
                auto su = cross(direction, u - p), sv = cross(direction, v - p);
                if (su == 0)
                    expected.push_back(rational(u));
                if (su * sv < 0)
                    expected.push_back(rational(u) + rational(v - u) * R(su, su - sv));
            }
            assert(intersection.has_value() == !expected.empty());
            if (intersection) {
                auto key = [&](const RP &v) { return dot(rational(direction), v); };
                auto [first, last] = std::minmax_element(expected.begin(), expected.end(),
                    [&](const RP &u, const RP &v) { return key(u) < key(v); });
                assert(intersection->s == *first && intersection->t == *last);
            }
        }

        auto width = convex_minimum_width(a);
        auto rectangle = convex_minimum_area_rectangle(a);
        assert(width.has_value() == !a.empty() && rectangle.has_value() == !a.empty());
        if (a.empty())
            continue;
        long double expected_width = 0, expected_area = 0;
        for (std::size_t i = 0; a.size() >= 2 && i < a.size(); ++i) {
            auto edge = a[(i + 1) % a.size()] - a[i];
            Int top = 0, left = 0, right = 0;
            for (auto p : a) {
                auto diff = p - a[i];
                top = std::max(top, cross(edge, diff));
                left = std::min(left, dot(edge, diff));
                right = std::max(right, dot(edge, diff));
            }
            auto w = static_cast<long double>(top) / std::sqrt(static_cast<long double>(norm(edge)));
            auto s = static_cast<long double>(top) * (right - left) / norm(edge);
            expected_width = i == 0 ? w : std::min(expected_width, w);
            expected_area = i == 0 ? s : std::min(expected_area, s);
        }
        assert(std::abs(*width - expected_width) < 1e-10L);
        assert(std::abs(rectangle->area - expected_area) < 1e-10L);
        if (rectangle->area > 0) {
            for (auto p : a)
                for (std::size_t i = 0; i < 4; ++i) {
                    auto u = rectangle->corners[i], v = rectangle->corners[(i + 1) % 4];
                    assert(cross(v - u, Point<long double>{static_cast<long double>(p.x), static_cast<long double>(p.y)} - u)
                           >= -1e-10L);
                }
        }
    }

    auto large = ConvexPolygon<double>{{0, 0}, {2e100, 0}, {2e100, 1e100}, {0, 1e100}};
    assert(std::abs(*convex_minimum_width(large) / 1e100 - 1) < 1e-14);
    assert(std::abs(convex_minimum_area_rectangle(large)->area / 2e200 - 1) < 1e-14);

    auto lattice = lattice_polygon_points(Polygon<Int>{{0, 0}, {0, 4}, {4, 4}, {4, 0}});
    assert(lattice.interior == 9 && lattice.boundary == 16);
    Int lo = std::numeric_limits<Int>::min(), hi = std::numeric_limits<Int>::max();
    auto huge = lattice_polygon_points(Polygon<Int>{{lo, lo}, {hi, lo}, {hi, hi}, {lo, hi}});
    BigInt span = BigInt(hi) - lo;
    assert(huge.interior == (span - 1) * (span - 1) && huge.boundary == span * 4);
}

void test_range_queries(std::mt19937 &random) {
    for (int round = 0; round < 120; ++round) {
        auto vertices = random_points(random, 8, 4), points = random_points(random, random() % 140, 6);
        TrianglePointCounter counter(vertices, points);
        PartitionTree<Int> unweighted(points);
        for (std::size_t i = 0; i < vertices.size(); ++i)
            for (std::size_t j = 0; j < vertices.size(); ++j)
                for (std::size_t k = 0; k < vertices.size(); ++k) {
                    auto a = vertices[i], b = vertices[j], c = vertices[k];
                    auto poly = hull({a, b, c});
                    Int strict = 0, closed = 0;
                    for (P p : points) {
                        strict += poly.contains(p, true);
                        closed += poly.contains(p);
                    }
                    assert(counter.count(i, j, k) == strict);
                    assert(counter.count(i, j, k, false) == closed);
                    assert(unweighted.count_triangle(a, b, c) == std::size_t(strict));
                    assert(unweighted.count_triangle(a, b, c, false) == std::size_t(closed));
                }

        std::vector<Int> weights;
        for (std::size_t i = 0; i < points.size(); ++i)
            weights.push_back(random_coordinate(random, 500));
        PartitionTree<Int> tree(points, weights);
        for (int q = 0; q < 30; ++q) {
            Int a = random_coordinate(random), b = random_coordinate(random), c = random_coordinate(random, 150);
            std::vector<std::size_t> expected;
            Int sum = 0;
            for (std::size_t i = 0; i < points.size(); ++i)
                if (a * points[i].x + b * points[i].y <= c) {
                    expected.push_back(i);
                    sum += weights[i];
                }
            auto ids = tree.report(a, b, c);
            std::sort(ids.begin(), ids.end());
            assert(ids == expected && tree.count(a, b, c) == ids.size() && tree.sum(a, b, c) == sum);
        }
    }

    const Int lo = std::numeric_limits<Int>::min(), hi = std::numeric_limits<Int>::max();
    const std::vector<P> corners{{lo, lo}, {hi, lo}, {lo, hi}, {hi, hi}};
    const std::vector<P> points{{0, -1}, {0, 0}, {-1, 0}, {-1, -1}, {0, 1}, {0, 1}, {hi, lo}};
    TrianglePointCounter counter(corners, points);
    assert(counter.count(0, 1, 2) == 1 && counter.count(0, 1, 2, false) == 4);
    PartitionTree<Int> tree(points);
    assert(tree.count_triangle(corners[0], corners[1], corners[2]) == 1);
    assert(tree.count_triangle(corners[0], corners[1], corners[2], false) == 4);
}

void test_dynamic_convex_hull(std::mt19937 &random) {
    DynamicConvexHull<Int> dynamic;
    std::vector<P> active;
    for (int i = 0; i < 14000; ++i) {
        if (active.empty() || random() % 100 < 51) {
            P p{random_coordinate(random, 6), random_coordinate(random, 6)};
            active.push_back(p);
            dynamic.insert(p);
        } else {
            auto next = active.begin() + random() % active.size();
            assert(dynamic.erase(*next));
            active.erase(next);
        }
        assert(!dynamic.erase({50, 50}) && dynamic.size() == active.size());
        assert(dynamic.hull() == hull(active));
        for (auto direction : random_points(random, 4, 5)) {
            const auto point = dynamic.support_point(direction);
            assert(point.has_value() == !active.empty());
            if (point) {
                auto expected = std::max_element(active.begin(), active.end(),
                    [&](P a, P b) { return dot(a, direction) < dot(b, direction); });
                assert(dot(*point, direction) == dot(*expected, direction));
            }
        }
    }
}

void test_delaunay_and_voronoi(std::mt19937 &random) {
    auto circle = voronoi_diagram(std::vector<P>{{1, 0}, {0, 1}, {-1, 0}, {0, -1}, {-1, 0}});
    assert(circle.vertices.size() == 1 && circle.edges.size() == 4);
    assert(circle.representative[4] == 2 && circle.cell(4) == circle.cell(2));
    assert(VoronoiNearestNeighbor(circle).nearest(P{0, 0}).has_value());
    for (Int i = 0; i < Int(circle.edges.size()); ++i) {
        const auto [a, b] = circle.edges[i].sites;
        auto anchor = circle.edge_anchor(i);
        assert(norm(anchor - circle.sites[a]) == norm(anchor - circle.sites[b]));
        assert(dot(circle.edge_direction(i), circle.sites[b] - circle.sites[a]) == R(0));
    }
    auto lo = std::numeric_limits<Int>::min(), hi = std::numeric_limits<Int>::max();
    VoronoiNearestNeighbor extreme(std::vector<P>{{lo, lo}, {hi, lo}, {lo, hi}, {hi, hi}});
    assert(extreme.nearest(P{0, 0}) == 3);
    assert(extreme.nearest(P{-1, -1}) == 0);
    const std::vector<Point<double>> floating{{0.25, 0.125}, {-0.25, 0.125}, {0, -0.125}, {0, -0.125}};
    VoronoiNearestNeighbor floating_search(floating);
    assert(floating_search.nearest(Point<double>{1e250, 0.125}) == 0);
    assert(floating_search.nearest(Point<double>{-1e250, 0.125}) == 1);
    assert(floating_search.nearest(Point<double>{0, -1e250}) == 2);

    for (int round = 0; round < 220; ++round) {
        auto points = random_points(random, random() % 47);
        for (auto &p : points) {
            if (round % 9 == 0)
                p.y = 0;
            if (round % 9 == 1)
                p.x = 0;
            if (round % 9 == 2)
                p.x = p.y;
        }
        auto mesh = delaunay_triangulation(points);
        std::set<std::array<Int, 2>> seen;
        for (std::size_t i = 0; i < mesh.edges.size(); ++i) {
            auto [u, v] = mesh.edges[i];
            assert(seen.insert({std::min(u, v), std::max(u, v)}).second);
            auto [left, right] = mesh.edge_faces[i];
            auto check_face = [&](Int a, Int b, Int face, Int other) {
                if (face < 0)
                    return;
                auto triangle = mesh.triangles[face];
                auto j = std::find(triangle.begin(), triangle.end(), a) - triangle.begin();
                assert(j < 3 && triangle[(j + 1) % 3] == b && mesh.neighbors[face][j] == other);
            };
            check_face(u, v, left, right);
            check_face(v, u, right, left);
        }
        for (auto [a, b, c] : mesh.triangles) {
            assert(cross(points[b] - points[a], points[c] - points[a]) > 0);
            for (auto q : points) {
                auto u = points[a] - q, v = points[b] - q, w = points[c] - q;
                assert(norm(u) * cross(v, w) - norm(v) * cross(u, w) + norm(w) * cross(u, v) <= 0);
            }
        }

        VoronoiNearestNeighbor search(points);
        for (auto q : random_points(random, 90, 60)) {
            auto nearest = search.nearest(q);
            if (points.empty()) {
                assert(!nearest);
            } else {
                assert(nearest && mesh.representative[*nearest] == *nearest);
                assert(std::ranges::none_of(points, [&](auto p) { return norm(p - q) < norm(points[*nearest] - q); }));
            }
        }
    }
}
}

int main() {
    auto polygon = Polygon<Int>{{0, 0}, {6, 0}, {6, 6}, {4, 6}, {4, 2}, {2, 2}, {2, 6}, {0, 6}};
    PersistentPolygonPointLocator<Int> locator(polygon);
    for (Int x = -1; x <= 7; ++x)
        for (Int y = -1; y <= 7; ++y)
            assert(locator.contains(P{x, y}) == contains(polygon, P{x, y}));

    std::mt19937 random(945123);
    test_convex_operations(random);
    test_range_queries(random);
    test_dynamic_convex_hull(random);
    test_delaunay_and_voronoi(random);
    std::cout << "advanced geometry: all cases passed\n";
}
