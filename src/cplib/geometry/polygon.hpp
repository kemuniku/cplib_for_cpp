#pragma once
#include <cplib/geometry/ccw.hpp>
#include <cplib/math/fractions.hpp>

namespace cplib {
template <class T> struct Polygon {
    std::vector<Point<T>> v;

    // 頂点数を返す。計算量 O(1)。
    Int len() const {
        return v.size();
    }

    // 頂点数を返す。計算量 O(1)。
    Int size() const {
        return len();
    }

    // 頂点列の先頭イテレータを返す。計算量 O(1)。
    auto begin() const {
        return v.begin();
    }

    // 頂点列の終端イテレータを返す。計算量 O(1)。
    auto end() const {
        return v.end();
    }
};

// 多角形の頂点数を返す。計算量 O(1)。
template <class T> Int len(const Polygon<T> &p) {
    return p.len();
}

// N個の頂点列 v をコピーして多角形を初期化する。計算量 O(N)。
template <class T> auto initPolygon(std::span<const Point<T>> v) {
    return Polygon<T>{{v.begin(), v.end()}};
}

// N個の頂点列をコピーして多角形を初期化する。計算量 O(N)。
template <class T> auto initPolygon(const std::vector<Point<T>> &v) {
    return Polygon<T>{v};
}

// N個の頂点列をコピーして多角形を初期化する。計算量 O(N)。
template <class T> auto initPolygon(std::initializer_list<Point<T>> v) {
    return Polygon<T>{{v.begin(), v.end()}};
}

// 整数版は符号付き面積の2倍、それ以外は面積を返す。O(N)。
// 符号付き面積を返す。整数型では面積の2倍となり、頂点列が時計回りなら負になる。
template <class T> T area(const Polygon<T> &p) {
    T out = T(0);
    for (Int i = 0; i < p.len(); ++i) {
        T c = cross(p.v[i], p.v[(i + 1) % p.len()]);
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
    for (Int i = 0; i < p.len(); ++i) {
        Int c = ccw(p.v[i], p.v[(i + 1) % p.len()], p.v[(i + 2) % p.len()]);
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
    std::reverse(reverse.v.begin(), reverse.v.end());
    return is_convex_ccw(reverse, strict);
}

// 点 p が多角形 poly の辺上にあるかどうかを判定する。頂点数をNとして計算量 O(N)。
// 隣り合う頂点は異なること。
template <class T> bool on_edge(const Polygon<T> &poly, Point<T> p) {
    for (Int i = 0; i < poly.len(); ++i)
        if (ccw(poly.v[i], poly.v[(i + 1) % poly.len()], p) == ON_SEGMENT)
            return true;
    return false;
}

// 頂点数をNとして計算量 O(N)。自己交差のない多角形 poly 内に点 p が存在するかを判定、辺上にある場合を含めない場合は strict = true を設定。
// 隣り合う頂点は異なること。
template <class T> bool contains(const Polygon<T> &poly, Point<T> p, bool strict = false) {
    if (on_edge(poly, p))
        return !strict;
    bool out = false;
    for (Int i = 0; i < poly.len(); ++i) {
        auto a = poly.v[i] - p, b = poly.v[(i + 1) % poly.len()] - p;
        if (a.y > b.y)
            std::swap(a, b);
        if (geometry_le(a.y, 0) && geometry_gt(b.y, 0) && geometry_lt(cross(a, b), 0))
            out = !out;
    }
    return out;
}

// 単調鎖法で凸包を返す。元実装どおり小入力の順序と、最後に先頭要素を消す出力順を保持。O(N log N)。
// N>=3の場合は重複点を含まないこと。
template <class T> Polygon<T> convex_hull(std::span<const Point<T>> v, bool strict = true) {
    Int n = v.size();
    if (n < 3)
        return {{v.begin(), v.end()}};
    std::vector<Point<T>> s(v.begin(), v.end());
    std::stable_sort(s.begin(), s.end());
    std::vector<Point<T>> vi = {s[0], s[1]};
    for (Int i = 2; i < n; ++i) {
        while (vi.size() >= 2 &&
               (strict ? ccw(vi[vi.size() - 2], vi.back(), s[i]) != COUNTER_CLOCKWISE
                       : ccw(vi[vi.size() - 2], vi.back(), s[i]) == CLOCKWISE))
            vi.pop_back();
        vi.push_back(s[i]);
    }
    std::size_t lower = vi.size();
    for (Int i = n - 2; i >= 0; --i) {
        while (vi.size() > lower &&
               (strict ? ccw(vi[vi.size() - 2], vi.back(), s[i]) != COUNTER_CLOCKWISE
                       : ccw(vi[vi.size() - 2], vi.back(), s[i]) == CLOCKWISE))
            vi.pop_back();
        vi.push_back(s[i]);
    }
    vi.erase(vi.begin());
    return {std::move(vi)};
}

// N個の点の凸包を返す。計算量 O(N log N)。
template <class T> auto convex_hull(const std::vector<Point<T>> &v, bool strict = true) {
    return convex_hull<T>(std::span<const Point<T>>(v), strict);
}
}
