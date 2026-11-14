#pragma once
#include <cplib/geometry/polygon.hpp>

namespace cplib {
namespace detail {
template <class T> Int convexSign(T value) {
    return geometry_lt(value, T(0)) ? -1 : geometry_gt(value, T(0)) ? 1 : 0;
}
template <class T> bool convexOnSegment(Point<T> a, Point<T> b, Point<T> p) {
    if (exact_equal(a, b))
        return p == a;
    return convexSign(cross(b - a, p - a)) == 0 && geometry_le(dot(p - a, p - b), T(0));
}
}
template <class T> class ConvexPolygon : public std::vector<Point<T>> {
    using Base = std::vector<Point<T>>;
    Int location(Point<T> p) const {
        const Base &vertices = *this;
        Int n = this->size();
        if (n == 0)
            return -1;
        if (n == 1)
            return p == vertices[0] ? 0 : -1;
        if (n == 2)
            return detail::convexOnSegment(vertices[0], vertices[1], p) ? 0 : -1;
        auto origin = vertices[0];
        Int first = detail::convexSign(cross(vertices[1] - origin, p - origin));
        Int last = detail::convexSign(cross(vertices.back() - origin, p - origin));
        if (first < 0 || last > 0)
            return -1;
        if (first == 0)
            return detail::convexOnSegment(origin, vertices[1], p) ? 0 : -1;
        if (last == 0)
            return detail::convexOnSegment(origin, vertices.back(), p) ? 0 : -1;
        Int lo = 1, hi = n - 1;
        while (hi - lo > 1) {
            Int mid = (lo + hi) / 2;
            if (cross(vertices[mid] - origin, p - origin) >= T(0))
                lo = mid;
            else
                hi = mid;
        }
        Int side = detail::convexSign(cross(vertices[hi] - vertices[lo], p - vertices[lo]));
        if (side == 0 && !detail::convexOnSegment(vertices[lo], vertices[hi], p))
            return -1;
        return side;
    }
public:
    using Base::Base;
    using Base::operator=;

    ConvexPolygon() = default;

    // 凸多角形の周回順頂点列を正規化する。共線・重複点も許す。時間・領域 O(N)。
    explicit ConvexPolygon(const std::vector<Point<T>> &input) {
        Base &vertices = *this;
        std::vector<Point<T>> points;
        for (auto p : input)
            if (points.empty() || !exact_equal(points.back(), p))
                points.push_back(p);
        if (points.size() > 1 && exact_equal(points.front(), points.back()))
            points.pop_back();
        if (points.size() <= 1) {
            vertices = std::move(points);
            return;
        }
        Int lo = std::min_element(points.begin(), points.end()) - points.begin();
        Int hi = std::max_element(points.begin(), points.end()) - points.begin();
        T orientation = T(0);
        for (Int i = 0; i < Int(points.size()); ++i) {
            T turn = cross(points[i] - points[lo], points[(i + 1) % points.size()] - points[lo]);
            if ((turn < T(0) ? -turn : turn) > (orientation < T(0) ? -orientation : orientation))
                orientation = turn;
        }
        if (orientation == T(0)) {
            vertices = {points[lo], points[hi]};
            return;
        }
        if (orientation < T(0)) {
            std::reverse(points.begin(), points.end());
            lo = Int(points.size()) - 1 - lo;
        }
        for (Int i = 0; i < Int(points.size()); ++i) {
            auto p = points[(lo + i) % points.size()];
            while (vertices.size() >= 2 &&
                   cross(vertices.back() - vertices[vertices.size() - 2],
                         p - vertices[vertices.size() - 2]) <= T(0))
                vertices.pop_back();
            vertices.push_back(p);
        }
        while (vertices.size() >= 3 &&
               cross(vertices.back() - vertices[vertices.size() - 2],
                     vertices.front() - vertices[vertices.size() - 2]) <= T(0))
            vertices.pop_back();
    }
    // Polygonをコピーして正規化する。時間・領域 O(N)。
    explicit ConvexPolygon(const Polygon<T> &poly) : ConvexPolygon(static_cast<const std::vector<Point<T>> &>(poly)) {
    }
    // 独立したPolygonを返す。時間・領域 O(N)。
    Polygon<T> toPolygon() const {
        return Polygon<T>(static_cast<const Base &>(*this));
    }
    // 点の包含を判定する。strict=trueは境界を除く。時間 O(log(N+1))。
    bool contains(Point<T> p, bool strict = false) const {
        return strict ? location(p) > 0 : location(p) >= 0;
    }
    // 境界上かを判定する。時間 O(log(N+1))。
    bool on_edge(Point<T> p) const {
        return location(p) == 0;
    }
    // 二乗直径を返す。空の多角形は不可。回転calipersで時間 O(N)、追加領域 O(1)。
    T diameter() const {
        const Base &vertices = *this;
        Int n = this->size();
        assert(n > 0);
        T best = T(0);
        auto consider = [&](Int i, Int j) {
            T distance = norm(vertices[i] - vertices[j]);
            if (distance > best)
                best = distance;
        };
        if (n == 2)
            consider(0, 1);
        if (n >= 3) {
            Int j = 1;
            for (Int i = 0; i < n; ++i) {
                Int next = (i + 1) % n;
                auto edge = vertices[next] - vertices[i];
                while (cross(edge, vertices[(j + 1) % n] - vertices[j]) > T(0))
                    j = (j + 1) % n;
                consider(i, j);
                consider(next, j);
                if (cross(edge, vertices[(j + 1) % n] - vertices[j]) == T(0)) {
                    consider(i, (j + 1) % n);
                    consider(next, (j + 1) % n);
                }
            }
        }
        return best;
    }
};
// 凸頂点列を前処理する。時間・領域 O(N)。
template <class T> auto initConvexPolygon(const std::vector<Point<T>> &v) {
    return ConvexPolygon<T>(v);
}
// 凸Polygonを前処理する。時間・領域 O(N)。
template <class T> auto initConvexPolygon(const Polygon<T> &poly) {
    return ConvexPolygon<T>(poly);
}
// 凸多角形の包含判定。時間 O(log(N+1))。
template <class T> bool contains(const ConvexPolygon<T> &poly, Point<T> p, bool strict = false) {
    return poly.contains(p, strict);
}
// 凸多角形の境界判定。時間 O(log(N+1))。
template <class T> bool on_edge(const ConvexPolygon<T> &poly, Point<T> p) {
    return poly.on_edge(p);
}
// 二乗直径を返す。空の多角形は不可。時間 O(N)、追加領域 O(1)。
template <class T> T diameter(const ConvexPolygon<T> &poly) {
    return poly.diameter();
}
// 凸Polygonの二乗直径を返す。空の多角形は不可。時間・領域 O(N)。
template <class T> T diameter(const Polygon<T> &poly) {
    return initConvexPolygon(poly).diameter();
}
// 正規化した頂点列をコピーする。時間・領域 O(N)。
template <class T> auto toPolygon(const ConvexPolygon<T> &poly) {
    return poly.toPolygon();
}
// 有向直線の左側と境界を残す。浮動小数点またはFractionを使う。時間・出力領域 O(N)。
template <class T> requires (!std::is_integral_v<T>)
Polygon<T> convex_cut(const Polygon<T> &poly, Line<T> line) {
    assert(line.s != line.t);
    Polygon<T> result;
    auto append = [&](Point<T> p) {
        if (result.empty() || result.back() != p)
            result.push_back(p);
    };
    const Int n = static_cast<Int>(poly.size());
    for (Int i = 0; i < n; ++i) {
        auto a = poly[i], b = poly[(i + 1) % n];
        T ca = cross(toVector(line), a - line.s), cb = cross(toVector(line), b - line.s);
        Int sa = detail::convexSign(ca), sb = detail::convexSign(cb);
        if (sa >= 0)
            append(a);
        if (sa * sb < 0)
            append(a + (b - a) * (ca / (ca - cb)));
    }
    if (result.size() > 1 && result.front() == result.back())
        result.pop_back();
    return result;
}
// 前処理済み凸多角形を切断する。時間・領域 O(N)。
template <class T> requires (!std::is_integral_v<T>)
auto convex_cut(const ConvexPolygon<T> &poly, Line<T> line) {
    return convex_cut(poly.toPolygon(), line);
}
}
