#pragma once
#include <cplib/geometry/base.hpp>
#include <cplib/math/bigint.hpp>
#include <bit>
#include <stdexcept>
#include <tuple>

namespace cplib {

struct DelaunayTriangulation {
    // すべての頂点番号は元の入力列の添字。同座標では最初の添字を使う。
    std::vector<Int> representative;
    std::vector<std::array<Int, 2>> edges;
    std::vector<std::array<Int, 3>> triangles;  // 反時計回り
    std::vector<std::array<Int, 2>> edge_faces; // edges[i] の有向辺に対する左、右の面。外部は -1。
    std::vector<std::array<Int, 3>> neighbors; // triangles[i][j] から次の頂点への辺を挟む面。
};

namespace detail {
template <class Coord> struct DelaunaySite {
    Coord x, y;
    Int id;
};

template <class Coord> struct DelaunayMesh {
    std::vector<DelaunaySite<Coord>> sites;
    DelaunayTriangulation topology; // 頂点番号は sites の添字。
    std::vector<std::pair<Int, Int>> duplicates;
};

template <bool WithFaces, class Coord>
DelaunayMesh<Coord> buildDelaunayMesh(std::vector<DelaunaySite<Coord>> sorted) {
    using Determinant = std::conditional_t<std::is_same_v<Coord, Int>, __int128, Coord>;
    DelaunayMesh<Coord> result;
    auto &sites = result.sites;
    auto &topology = result.topology;
    std::sort(sorted.begin(), sorted.end(),
              [](const auto &a, const auto &b) { return std::tie(a.x, a.y, a.id) < std::tie(b.x, b.y, b.id); });
    sites.reserve(sorted.size());
    topology.representative.resize(sorted.size());
    for (auto &p : sorted) {
        if (!sites.empty() && sites.back().x == p.x && sites.back().y == p.y)
            result.duplicates.emplace_back(sites.back().id, p.id);
        else
            sites.push_back(std::move(p));
        topology.representative[p.id] = sites.size() - 1;
    }
    if (sites.size() <= 1)
        return result;

    struct Quad {
        std::array<Int, 4> next;
        std::array<Int, 2> vertex;
    };

    std::vector<Quad> edges;
    edges.reserve(3 * sites.size());
    std::vector<Int> unused;
    auto rotate = [](Int e) { return (e & ~Int(3)) | ((e + 1) & 3); };
    auto onext = [&](Int e) { return edges[e >> 2].next[e & 3]; };
    auto origin = [&](Int e) { return edges[e >> 2].vertex[(e >> 1) & 1]; };
    auto destination = [&](Int e) { return origin(e ^ 2); };
    auto oprev = [&](Int e) { return rotate(onext(rotate(e))); };
    auto lnext = [&](Int e) { return rotate(onext(rotate(e) ^ 2)); };
    auto splice = [&](Int a, Int b) {
        Int alpha = rotate(onext(a)), beta = rotate(onext(b));
        std::swap(edges[a >> 2].next[a & 3], edges[b >> 2].next[b & 3]);
        std::swap(edges[alpha >> 2].next[alpha & 3], edges[beta >> 2].next[beta & 3]);
    };
    auto makeEdge = [&](Int u, Int v) {
        Int index;
        if (unused.empty()) {
            index = edges.size();
            edges.emplace_back();
        } else {
            index = unused.back();
            unused.pop_back();
        }
        Int e = index << 2;
        edges[index] = {{{e, e + 3, e + 2, e + 1}}, {{u, v}}};
        return e;
    };
    auto connect = [&](Int a, Int b) {
        Int e = makeEdge(destination(a), origin(b));
        splice(e, lnext(a));
        splice(e ^ 2, b);
        return e;
    };
    auto remove = [&](Int e) {
        splice(e, oprev(e));
        splice(e ^ 2, oprev(e ^ 2));
        edges[e >> 2].vertex[0] = -1;
        unused.push_back(e >> 2);
    };
    auto orientation = [&](Int a, Int b, Int c) {
        return (sites[b].x - sites[a].x) * (sites[c].y - sites[a].y) -
               (sites[b].y - sites[a].y) * (sites[c].x - sites[a].x);
    };
    auto leftOf = [&](Int p, Int e) { return orientation(origin(e), destination(e), p) > 0; };
    auto rightOf = [&](Int p, Int e) { return orientation(origin(e), destination(e), p) < 0; };
    auto inCircle = [&](Int a, Int b, Int c, Int d) {
        Coord ax = sites[a].x - sites[d].x, ay = sites[a].y - sites[d].y,
              bx = sites[b].x - sites[d].x, by = sites[b].y - sites[d].y,
              cx = sites[c].x - sites[d].x, cy = sites[c].y - sites[d].y;
        Determinant det = Determinant(ax * ax + ay * ay) * (bx * cy - by * cx) -
                          Determinant(bx * bx + by * by) * (ax * cy - ay * cx) +
                          Determinant(cx * cx + cy * cy) * (ax * by - ay * bx);
        return det > 0;
    };
    auto triangulate = [&](auto &&self, Int first, Int last) -> std::pair<Int, Int> {
        if (last - first <= 3) {
            Int a = makeEdge(first, first + 1);
            if (last - first == 2)
                return {a, a ^ 2};
            Int b = makeEdge(first + 1, first + 2);
            splice(a ^ 2, b);
            Coord turn = orientation(first, first + 1, first + 2);
            if (turn == 0)
                return {a, b ^ 2};
            Int c = connect(b, a);
            return turn > 0 ? std::pair{a, b ^ 2} : std::pair{c ^ 2, c};
        }
        Int mid = (first + last) / 2;
        auto [lo, li] = self(self, first, mid);
        auto [ri, ro] = self(self, mid, last);
        while (true) {
            if (leftOf(origin(ri), li))
                li = lnext(li);
            else if (rightOf(origin(li), ri))
                ri = onext(ri ^ 2);
            else
                break;
        }
        Int base = connect(ri ^ 2, li);
        if (origin(li) == origin(lo))
            lo = base ^ 2;
        if (origin(ri) == origin(ro))
            ro = base;
        while (true) {
            Int left = onext(base ^ 2);
            if (rightOf(destination(left), base))
                while (inCircle(destination(base), origin(base), destination(left),
                                destination(onext(left)))) {
                    Int next = onext(left);
                    remove(left);
                    left = next;
                }
            Int right = oprev(base);
            if (rightOf(destination(right), base))
                while (inCircle(destination(base), origin(base), destination(right),
                                destination(oprev(right)))) {
                    Int next = oprev(right);
                    remove(right);
                    right = next;
                }
            bool lv = rightOf(destination(left), base), rv = rightOf(destination(right), base);
            if (!lv && !rv)
                break;
            if (!lv || (rv && inCircle(destination(left), origin(left), origin(right),
                                       destination(right))))
                base = connect(right, base ^ 2);
            else
                base = connect(base ^ 2, left ^ 2);
        }
        return {lo, ro};
    };
    triangulate(triangulate, 0, sites.size());

    if constexpr (WithFaces) {
        constexpr Int unvisited = -2;
        std::vector<Int> face_of(2 * edges.size(), unvisited);
        std::vector<Int> cycle;
        auto face_at = [&](Int e) -> Int & { return face_of[e >> 1]; };
        for (Int i = 0; i < Int(edges.size()); ++i) {
            if (edges[i].vertex[0] < 0)
                continue;
            for (Int start : {4 * i, 4 * i + 2}) {
                if (face_at(start) != unvisited)
                    continue;
                cycle.clear();
                Int e = start;
                do {
                    cycle.push_back(e);
                    e = lnext(e);
                } while (e != start);
                Int face = -1;
                if (cycle.size() == 3) {
                    auto a = origin(cycle[0]), b = origin(cycle[1]), c = origin(cycle[2]);
                    if (orientation(a, b, c) > 0) {
                        face = topology.triangles.size();
                        topology.triangles.push_back({a, b, c});
                    }
                }
                for (Int boundary : cycle)
                    face_at(boundary) = face;
            }
        }
        topology.neighbors.resize(topology.triangles.size());
        for (Int i = 0; i < Int(edges.size()); ++i) {
            if (edges[i].vertex[0] < 0)
                continue;
            topology.edges.push_back(edges[i].vertex);
            topology.edge_faces.push_back({face_at(4 * i), face_at(4 * i + 2)});
        }
        for (Int i = 0; i < Int(edges.size()); ++i) {
            if (edges[i].vertex[0] < 0)
                continue;
            for (Int e : {4 * i, 4 * i + 2}) {
                Int face = face_at(e);
                if (face < 0)
                    continue;
                const auto &triangle = topology.triangles[face];
                auto vertex = std::find(triangle.begin(), triangle.end(), origin(e));
                topology.neighbors[face][vertex - triangle.begin()] = face_at(e ^ 2);
            }
        }
    } else {
        for (const auto &e : edges)
            if (e.vertex[0] >= 0)
                topology.edges.push_back(e.vertex);
    }
    return result;
}

struct DelaunayBinaryCoordinate {
    __uint128_t significand = 0;
    int exponent = 0, upper_exponent = 0;
    bool negative = false;
};

template <FloatingPoint T> DelaunayBinaryCoordinate delaunayBinaryCoordinate(T value) {
    constexpr int digits = NumericLimits<T>::digits;
    static_assert(NumericLimits<T>::radix == 2 && 0 < digits && digits <= 128);
    if (!isfinite(value))
        throw std::invalid_argument("geometry requires finite coordinates");
    if (value == 0)
        return {};
    bool negative = value < 0;
    int exponent;
    T fraction = frexp(negative ? -value : value, &exponent);
    auto significand = static_cast<__uint128_t>(ldexp(fraction, digits));
    UInt low = static_cast<UInt>(significand);
    int trailing = low ? std::countr_zero(low)
                       : 64 + std::countr_zero(static_cast<UInt>(significand >> 64));
    return {significand >> trailing, exponent - digits + trailing, exponent, negative};
}

inline BigInt delaunaySignificand(DelaunayBinaryCoordinate value) {
    BigInt magnitude(static_cast<UInt>(value.significand));
    UInt high = static_cast<UInt>(value.significand >> 64);
    if (high != 0)
        magnitude += BigInt(high) << 64;
    return value.negative ? -magnitude : magnitude;
}

template <FloatingPoint T, class Id, class F>
auto withDelaunayBinaryCoordinates(std::span<const Point<T, Id>> points, F &&run) {
    int minExponent = std::numeric_limits<int>::max();
    int maxExponent = std::numeric_limits<int>::min();
    auto decompose = [&](T value) {
        auto coordinate = delaunayBinaryCoordinate(value);
        if (coordinate.significand != 0) {
            minExponent = std::min(minExponent, coordinate.exponent);
            maxExponent = std::max(maxExponent, coordinate.upper_exponent);
        }
        return coordinate;
    };

    std::vector<std::array<DelaunayBinaryCoordinate, 2>> coordinates;
    coordinates.reserve(points.size());
    for (const auto &point : points)
        coordinates.push_back({decompose(point.x), decompose(point.y)});

    // 共通の2のべきを掛けて整数にする。向き、外接円判定、辺の長さの順序は保たれる。
    auto apply = [&](auto scale) {
        using Coord = decltype(scale(DelaunayBinaryCoordinate{}));
        std::vector<DelaunaySite<Coord>> sites;
        sites.reserve(points.size());
        for (Int i = 0; i < Int(coordinates.size()); ++i) {
            const auto &[x, y] = coordinates[i];
            sites.push_back({scale(x), scale(y), i});
        }
        return run(std::move(sites));
    };
    // スケール後の絶対値が29ビット以下なら、整数版と同じ計算幅を使える。
    if (maxExponent == std::numeric_limits<int>::min() || maxExponent - minExponent <= 29) {
        return apply([&](const DelaunayBinaryCoordinate &value) -> Int {
            if (value.significand == 0)
                return 0;
            Int magnitude = static_cast<Int>(value.significand << (value.exponent - minExponent));
            return value.negative ? -magnitude : magnitude;
        });
    }
    return apply([&](const DelaunayBinaryCoordinate &value) -> BigInt {
        if (value.significand == 0)
            return {};
        return delaunaySignificand(value) << (value.exponent - minExponent);
    });
}

template <class T, class Id, class F>
auto withDelaunayCoordinates(std::span<const Point<T, Id>> points, F &&run) {
    if constexpr (FloatingPoint<T>) {
        return withDelaunayBinaryCoordinates(points, std::forward<F>(run));
    } else {
        static_assert(std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int));
        auto coordinate = [](T x) { return -1000000000 <= x && x <= 1000000000; };
        bool narrow = std::all_of(points.begin(), points.end(),
                                  [&](const auto &p) { return coordinate(p.x) && coordinate(p.y); });
        auto apply = [&]<class Coord>() {
            std::vector<DelaunaySite<Coord>> sites;
            sites.reserve(points.size());
            for (Int i = 0; i < Int(points.size()); ++i)
                sites.push_back({Coord(points[i].x), Coord(points[i].y), i});
            return run(std::move(sites));
        };
        if (narrow)
            return apply.template operator()<Int>();
        return apply.template operator()<BigInt>();
    }
}
} // namespace detail

// O(N log N) 回の算術・比較と O(N) 領域。全点共線なら隣り合う点の辺のみを返す。
// 符号付き整数または有限な2進浮動小数点に対応。必要に応じて BigInt で判定する。
template <class T, class Id = void>
DelaunayTriangulation delaunay_triangulation(std::span<const Point<T, Id>> points) {
    return detail::withDelaunayCoordinates(points, [](auto sites) {
        auto mesh = detail::buildDelaunayMesh<true>(std::move(sites));
        auto original = [&](Int i) { return mesh.sites[i].id; };
        for (Int &vertex : mesh.topology.representative)
            vertex = original(vertex);
        for (auto &edge : mesh.topology.edges)
            for (Int &vertex : edge)
                vertex = original(vertex);
        for (auto &triangle : mesh.topology.triangles)
            for (Int &vertex : triangle)
                vertex = original(vertex);
        return std::move(mesh.topology);
    });
}

template <class T, class Id = void>
DelaunayTriangulation delaunay_triangulation(const std::vector<Point<T, Id>> &points) {
    return delaunay_triangulation(std::span<const Point<T, Id>>(points));
}

}
