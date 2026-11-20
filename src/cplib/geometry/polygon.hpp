#pragma once
#include <cplib/geometry/ccw.hpp>

namespace cplib {
template <class T, class Id = void> struct Polygon : std::vector<Point<T, Id>> {
    using Base = std::vector<Point<T, Id>>;
    using Base::Base;
    using Base::operator=;

    Polygon() = default;
    Polygon(const Base &v) : Base(v) {}
    Polygon(Base &&v) : Base(std::move(v)) {}
};

// N個の頂点列 v をコピーして多角形を初期化する。計算量 O(N)。
// 整数幅: 頂点のコピーのみ。B=M、M=10^9ならint可。
template <class T, class Id = void> auto initPolygon(std::span<const Point<T, Id>> v) {
    return Polygon<T, Id>(v.begin(), v.end());
}

// N個の頂点列をコピーして多角形を初期化する。計算量 O(N)。
// 整数幅: 頂点のコピーのみ。B=M、M=10^9ならint可。
template <class T, class Id = void> auto initPolygon(const std::vector<Point<T, Id>> &v) {
    return Polygon<T, Id>(v);
}

// N個の頂点列をコピーして多角形を初期化する。計算量 O(N)。
// 整数幅: 頂点のコピーのみ。B=M、M=10^9ならint可。
template <class T, class Id = void> auto initPolygon(std::initializer_list<Point<T, Id>> v) {
    return Polygon<T, Id>(v);
}

// 整数版は符号付き面積の2倍、それ以外は面積を返す。O(N)。
// 符号付き面積を返す。整数型では面積の2倍となり、頂点列が時計回りなら負になる。
// 整数幅: B=2*max(1,N)*M^2。累積和も含めW(B)ビット（Nに応じて追加幅が必要）。
// M=10^9ではN<=4ならlong longで保証、一般のNでは__int128等を選ぶ。
// 最終的な面積だけでなく、各辺の外積を加算する途中の値も収まること。
template <class T, class Id = void> T area(const Polygon<T, Id> &p) {
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
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> bool is_convex_ccw(const Polygon<T, Id> &p, bool strict = true) {
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
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> bool is_convex(const Polygon<T, Id> &p, bool strict = true) {
    if (is_convex_ccw(p, strict))
        return true;
    auto reverse = p;
    std::reverse(reverse.begin(), reverse.end());
    return is_convex_ccw(reverse, strict);
}

// 点 p が多角形 poly の辺上にあるかどうかを判定する。頂点数をNとして計算量 O(N)。
// 隣り合う頂点は異なること。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> bool on_edge(const Polygon<T, Id> &poly, Point<T> p) {
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
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> Int contains(const Polygon<T, Id> &poly, Point<T> p) {
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
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> bool contains(const Polygon<T, Id> &poly, Point<T> p, bool strict) {
    Int result = contains(poly, p);
    return strict ? result == 2 : result != 0;
}

namespace detail {
// 辞書順に並んだ点列から凸包を構築する。時間・領域 O(N)。
template <class T, class Id> Polygon<T, Id> convexHullFromSorted(std::vector<Point<T, Id>> points, bool strict) {
    points.erase(std::unique(points.begin(), points.end(), [](const auto &a, const auto &b) {
        return exact_equal(a, b);
    }), points.end());
    if (points.size() <= 2)
        return Polygon<T, Id>(std::move(points));
    const T zero = T(0);
    bool collinear = true;
    for (std::size_t i = 1; i + 1 < points.size(); ++i)
        if (cross(points.back() - points.front(), points[i] - points.front()) != zero)
            collinear = false;
    if (collinear) {
        if (strict)
            return Polygon<T, Id>{points.front(), points.back()};
        return Polygon<T, Id>(std::move(points));
    }
    std::vector<Point<T, Id>> lower, upper;
    auto append = [&](auto &chain, const Point<T, Id> &p) {
        while (chain.size() >= 2) {
            T turn = cross(chain.back() - chain[chain.size() - 2], p - chain[chain.size() - 2]);
            if (strict ? turn > zero : turn >= zero)
                break;
            chain.pop_back();
        }
        chain.push_back(p);
    };
    for (const auto &p : points)
        append(lower, p);
    for (auto it = points.rbegin(); it != points.rend(); ++it)
        append(upper, *it);
    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());
    return Polygon<T, Id>(std::move(lower));
}
}

// 辞書順最小点から反時計回りの凸包を返す。EPSを使わず重複点を除く。時間 O(N log N)、領域 O(N)。
// ID付きの同座標点は入力順で最初の点を残す。
// strict=falseは辺上の点も残す。全点共線ならstrict時は両端、それ以外は各点を一度返す。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> Polygon<T, Id> convex_hull(std::span<const Point<T, Id>> v, bool strict = true) {
    std::vector<Point<T, Id>> points(v.begin(), v.end());
    if constexpr (std::is_void_v<Id>)
        std::sort(points.begin(), points.end());
    else
        std::stable_sort(points.begin(), points.end());
    return detail::convexHullFromSorted(std::move(points), strict);
}

// N個の点の凸包を返す。計算量 O(N log N)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> auto convex_hull(const std::vector<Point<T, Id>> &v, bool strict = true) {
    return convex_hull<T, Id>(std::span<const Point<T, Id>>(v), strict);
}
}
