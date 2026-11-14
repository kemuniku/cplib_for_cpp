#pragma once
#include <cplib/geometry/ccw.hpp>

namespace cplib {
template <class T> struct Polygon : std::vector<Point<T>> {
    using Base = std::vector<Point<T>>;
    using Base::Base;
    using Base::operator=;

    Polygon() = default;
    Polygon(const Base &v) : Base(v) {}
    Polygon(Base &&v) : Base(std::move(v)) {}
};

// N個の頂点列 v をコピーして多角形を初期化する。計算量 O(N)。
template <class T> auto initPolygon(std::span<const Point<T>> v) {
    return Polygon<T>(v.begin(), v.end());
}

// N個の頂点列をコピーして多角形を初期化する。計算量 O(N)。
template <class T> auto initPolygon(const std::vector<Point<T>> &v) {
    return Polygon<T>(v);
}

// N個の頂点列をコピーして多角形を初期化する。計算量 O(N)。
template <class T> auto initPolygon(std::initializer_list<Point<T>> v) {
    return Polygon<T>(v);
}

// 整数版は符号付き面積の2倍、それ以外は面積を返す。O(N)。
// 符号付き面積を返す。整数型では面積の2倍となり、頂点列が時計回りなら負になる。
template <class T> T area(const Polygon<T> &p) {
    T out = T(0);
    const Int n = static_cast<Int>(p.size());
    for (Int i = 0; i < n; ++i) {
        T c = cross(p[i], p[(i + 1) % n]);
        if constexpr (std::is_integral_v<T>)
            out += c;
        else
            out += c / T(2);
    }
    return out;
}

// 頂点数をNとして計算量 O(N)。頂点列が反時計回りの多角形に対して、凸かどうかを判定する。
// 隣り合う頂点は異なること。
template <class T> bool is_convex_ccw(const Polygon<T> &p, bool strict = true) {
    const Int n = static_cast<Int>(p.size());
    for (Int i = 0; i < n; ++i) {
        Int c = ccw(p[i], p[(i + 1) % n], p[(i + 2) % n]);
        if (strict && c != COUNTER_CLOCKWISE)
            return false;
        if (!strict && c == CLOCKWISE)
            return false;
    }
    return true;
}

// 頂点数をNとして計算量 O(N)。多角形が凸かどうかを判定、ちょうど 180 度の内角を許容する場合は strict = false を設定。
// 隣り合う頂点は異なること。
template <class T> bool is_convex(const Polygon<T> &p, bool strict = true) {
    if (is_convex_ccw(p, strict))
        return true;
    auto reverse = p;
    std::reverse(reverse.begin(), reverse.end());
    return is_convex_ccw(reverse, strict);
}

// 点 p が多角形 poly の辺上にあるかどうかを判定する。頂点数をNとして計算量 O(N)。
// 隣り合う頂点は異なること。
template <class T> bool on_edge(const Polygon<T> &poly, Point<T> p) {
    const Int n = static_cast<Int>(poly.size());
    for (Int i = 0; i < n; ++i) {
        auto a = poly[i], b = poly[(i + 1) % n];
        if (exact_equal(a, b) || a == b) {
            if (p == a || p == b)
                return true;
        } else if (ccw(a, b, p) == ON_SEGMENT) {
            return true;
        }
    }
    return false;
}

// 頂点数をNとして計算量 O(N)。自己交差のない多角形 poly に対し、点 p が内部なら2、辺上なら1、外部なら0を返す。
// 隣り合う頂点は異なること。
template <class T> Int contains(const Polygon<T> &poly, Point<T> p) {
    if (on_edge(poly, p))
        return 1;
    bool out = false;
    const Int n = static_cast<Int>(poly.size());
    for (Int i = 0; i < n; ++i) {
        auto a = poly[i] - p, b = poly[(i + 1) % n] - p;
        if (a.y > b.y)
            std::swap(a, b);
        if (geometry_le(a.y, 0) && geometry_gt(b.y, 0) && geometry_lt(cross(a, b), 0))
            out = !out;
    }
    return out ? 2 : 0;
}

// strict=trueでは辺上を除く真偽値を返す。計算量 O(N)。
template <class T> bool contains(const Polygon<T> &poly, Point<T> p, bool strict) {
    Int result = contains(poly, p);
    return strict ? result == 2 : result != 0;
}

// 辞書順最小点から反時計回りの凸包を返す。EPSを使わず重複点を除く。時間 O(N log N)、領域 O(N)。
// strict=falseは辺上の点も残す。全点共線ならstrict時は両端、それ以外は各点を一度返す。
template <class T> Polygon<T> convex_hull(std::span<const Point<T>> v, bool strict = true) {
    std::vector<Point<T>> points(v.begin(), v.end());
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end(), [](auto a, auto b) {
        return exact_equal(a, b);
    }), points.end());
    if (points.size() <= 2)
        return Polygon<T>(std::move(points));
    const T zero = T(0);
    bool collinear = true;
    for (std::size_t i = 1; i + 1 < points.size(); ++i)
        if (cross(points.back() - points.front(), points[i] - points.front()) != zero)
            collinear = false;
    if (collinear) {
        if (strict)
            return Polygon<T>{points.front(), points.back()};
        return Polygon<T>(std::move(points));
    }
    std::vector<Point<T>> lower, upper;
    auto append = [&](auto &chain, Point<T> p) {
        while (chain.size() >= 2) {
            T turn = cross(chain.back() - chain[chain.size() - 2], p - chain[chain.size() - 2]);
            if (strict ? turn > zero : turn >= zero)
                break;
            chain.pop_back();
        }
        chain.push_back(p);
    };
    for (auto p : points)
        append(lower, p);
    for (auto it = points.rbegin(); it != points.rend(); ++it)
        append(upper, *it);
    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());
    return Polygon<T>(std::move(lower));
}

// N個の点の凸包を返す。計算量 O(N log N)。
template <class T> auto convex_hull(const std::vector<Point<T>> &v, bool strict = true) {
    return convex_hull<T>(std::span<const Point<T>>(v), strict);
}
}
