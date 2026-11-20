#pragma once
#include <cplib/geometry/base.hpp>
#include <cplib/collections/unionfind.hpp>
#include <cplib/math/bigint.hpp>
#include <bit>
#include <stdexcept>
#include <tuple>

namespace cplib {
namespace detail {
template <class Coord> struct EuclideanMstSite {
    Coord x, y;
    Int id;
};

template <class Coord>
std::vector<std::pair<Int, Int>> euclideanMst(std::vector<EuclideanMstSite<Coord>> sorted) {
    using Determinant = std::conditional_t<std::is_same_v<Coord, Int>, __int128, Coord>;
    std::vector<std::pair<Int, Int>> out;
    if (sorted.size() <= 1)
        return out;
    std::sort(sorted.begin(), sorted.end(),
              [](const auto &a, const auto &b) { return std::tie(a.x, a.y, a.id) < std::tie(b.x, b.y, b.id); });
    std::vector<EuclideanMstSite<Coord>> sites;
    sites.reserve(sorted.size());
    for (auto &p : sorted) {
        if (!sites.empty() && sites.back().x == p.x && sites.back().y == p.y)
            out.emplace_back(sites.back().id, p.id);
        else
            sites.push_back(std::move(p));
    }
    if (sites.size() == 1)
        return out;

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
    std::vector<std::tuple<Coord, Int, Int>> candidates;
    for (const auto &e : edges) {
        Int u = e.vertex[0], v = e.vertex[1];
        if (u < 0)
            continue;
        Coord dx = sites[u].x - sites[v].x, dy = sites[u].y - sites[v].y;
        candidates.emplace_back(dx * dx + dy * dy, u, v);
    }
    std::sort(candidates.begin(), candidates.end());
    UnionFind uf(sites.size());
    for (const auto &[weight, u, v] : candidates) {
        (void)weight;
        if (uf.issame(u, v))
            continue;
        uf.unite(u, v);
        out.emplace_back(std::min(sites[u].id, sites[v].id), std::max(sites[u].id, sites[v].id));
        if (uf.count == 1)
            break;
    }
    return out;
}

template <FloatingPoint T, class Id>
std::vector<std::pair<Int, Int>> euclideanMstFloating(std::span<const Point<T, Id>> points) {
    constexpr int digits = NumericLimits<T>::digits;
    static_assert(NumericLimits<T>::radix == 2 && 0 < digits && digits <= 128);

    struct BinaryCoordinate {
        __uint128_t significand = 0;
        int exponent = 0;
        bool negative = false;
    };

    int minExponent = std::numeric_limits<int>::max();
    int maxExponent = std::numeric_limits<int>::min();
    auto decompose = [&](T value) -> BinaryCoordinate {
        if (!isfinite(value))
            throw std::invalid_argument("euclidean_mst requires finite coordinates");
        if (value == 0)
            return {};
        bool negative = value < 0;
        int exponent;
        T fraction = frexp(negative ? -value : value, &exponent);
        auto significand = static_cast<__uint128_t>(ldexp(fraction, digits));
        UInt low = static_cast<UInt>(significand);
        int trailing = low ? std::countr_zero(low)
                           : 64 + std::countr_zero(static_cast<UInt>(significand >> 64));
        int normalizedExponent = exponent - digits + trailing;
        minExponent = std::min(minExponent, normalizedExponent);
        maxExponent = std::max(maxExponent, exponent);
        return {significand >> trailing, normalizedExponent, negative};
    };

    std::vector<std::array<BinaryCoordinate, 2>> coordinates;
    coordinates.reserve(points.size());
    for (const auto &point : points)
        coordinates.push_back({decompose(point.x), decompose(point.y)});

    // 共通の2のべきを掛けて整数にする。向き、外接円判定、辺の長さの順序は保たれる。
    auto run = [&](auto scale) {
        using Coord = decltype(scale(BinaryCoordinate{}));
        std::vector<EuclideanMstSite<Coord>> sites;
        sites.reserve(points.size());
        for (Int i = 0; i < Int(coordinates.size()); ++i) {
            const auto &[x, y] = coordinates[i];
            sites.push_back({scale(x), scale(y), i});
        }
        return euclideanMst(std::move(sites));
    };
    // スケール後の絶対値が29ビット以下なら、整数版と同じ計算幅を使える。
    if (maxExponent == std::numeric_limits<int>::min() || maxExponent - minExponent <= 29) {
        return run([&](const BinaryCoordinate &value) -> Int {
            if (value.significand == 0)
                return 0;
            Int magnitude = static_cast<Int>(value.significand << (value.exponent - minExponent));
            return value.negative ? -magnitude : magnitude;
        });
    }
    return run([&](const BinaryCoordinate &value) -> BigInt {
        if (value.significand == 0)
            return {};
        BigInt magnitude(static_cast<UInt>(value.significand));
        UInt high = static_cast<UInt>(value.significand >> 64);
        if (high != 0)
            magnitude += BigInt(high) << 64;
        if (value.exponent != minExponent)
            magnitude <<= value.exponent - minExponent;
        if (value.negative)
            return -magnitude;
        return magnitude;
    });
}
} // namespace detail

// 厳密なDelaunay三角形分割とKruskal法を使い、MSTを入力の頂点番号の組で返す。
// O(N log N)回の算術・比較、O(N)個の数値を格納。重複点も別頂点として扱い、N <= 1なら空列を返す。
// 符号付き整数座標は絶対値10^9以下。差・外積・二乗距離はInt（B=8M^2）、外接円判定は__int128（B=96M^4）。
// 有限な2進浮動小数点座標は表現された値を厳密に整数化する。多倍長演算の費用は指数の幅に依存する。
// 非有限の浮動小数点座標にはstd::invalid_argumentを送出する。
template <class T, class Id = void>
std::vector<std::pair<Int, Int>> euclidean_mst(std::span<const Point<T, Id>> points) {
    if constexpr (FloatingPoint<T>) {
        return detail::euclideanMstFloating(points);
    } else {
        static_assert(std::is_integral_v<T> && std::is_signed_v<T>);
        std::vector<detail::EuclideanMstSite<Int>> sites;
        sites.reserve(points.size());
        for (Int i = 0; i < Int(points.size()); ++i) {
            Int x = points[i].x, y = points[i].y;
            assert(-1000000000 <= x && x <= 1000000000 && -1000000000 <= y && y <= 1000000000);
            sites.push_back({x, y, i});
        }
        return detail::euclideanMst(std::move(sites));
    }
}

// 座標の制約はspan版と同じ。
template <class T, class Id = void> auto euclidean_mst(const std::vector<Point<T, Id>> &p) {
    return euclidean_mst<T, Id>(std::span<const Point<T, Id>>(p));
}

// 座標の組で表すN点のMSTを返す。座標の制約はPoint版と同じ。
template <class T> auto euclidean_mst(std::span<const std::pair<T, T>> p) {
    std::vector<Point<T>> v;
    v.reserve(p.size());
    for (auto [x, y] : p)
        v.push_back({x, y});
    return euclidean_mst(v);
}

// 座標の制約はPoint版と同じ。
template <class T> auto euclidean_mst(const std::vector<std::pair<T, T>> &p) {
    return euclidean_mst<T>(std::span<const std::pair<T, T>>(p));
}
}
