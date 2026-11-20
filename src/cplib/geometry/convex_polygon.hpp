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
template <class T, class Id = void> class ConvexPolygon : public std::vector<Point<T, Id>> {
    using Base = std::vector<Point<T, Id>>;
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
    // 入力条件: 凸多角形の頂点が時計回りまたは反時計回りの周回順に並んでいること。
    // 凸包構築・凸性検査は行わない。任意の点集合Pには initConvexPolygon(convex_hull(P)) を使う。
    // 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
    explicit ConvexPolygon(const std::vector<Point<T, Id>> &input) {
        Base &vertices = *this;
        Base points;
        for (const auto &p : input)
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
            const auto &p = points[(lo + i) % points.size()];
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
    // 入力条件: polyは凸多角形で、頂点が時計回りまたは反時計回りの周回順に並んでいること。
    // 凸包構築・凸性検査は行わない。任意の点集合Pには initConvexPolygon(convex_hull(P)) を使う。
    // 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
    explicit ConvexPolygon(const Polygon<T, Id> &poly) : ConvexPolygon(static_cast<const Base &>(poly)) {
    }
    // 独立したPolygonを返す。時間・領域 O(N)。
    // 整数幅: 座標のコピーのみ。B=M、M=10^9ならint可。
    Polygon<T, Id> toPolygon() const {
        return Polygon<T, Id>(static_cast<const Base &>(*this));
    }
    // 点の包含を判定する。strict=trueは境界を除く。時間 O(log(N+1))。
    // 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
    bool contains(Point<T> p, bool strict = false) const {
        return strict ? location(p) > 0 : location(p) >= 0;
    }
    // 境界上かを判定する。時間 O(log(N+1))。
    // 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
    bool on_edge(Point<T> p) const {
        return location(p) == 0;
    }
    // dot(direction, p) を最大にする頂点の、正規化後の添字を返す。同値ならいずれか1つ。
    // 正規化済みの空でない凸多角形が必要。方向がゼロなら0。時間 O(log(N+1))、追加領域 O(1)。
    // 整数幅: B=8M^2。Mは方向の成分の絶対値も含む。M=10^9ならT=long longで十分。
    Int support_index(Point<T> direction) const {
        const Base &vertices = *this;
        const Int n = static_cast<Int>(vertices.size());
        assert(n > 0);
        if (n == 1 || (direction.x == T(0) && direction.y == T(0)))
            return 0;
        if (n == 2)
            return dot(direction, vertices[1] - vertices[0]) > T(0) ? 1 : 0;

        // 最初の辺を基準に、内積が増加から非増加に変わる辺の方向を二分探索する。
        const Point<T> tangent{T(-direction.y), direction.x};
        const Point<T> anchor = vertices[1] - vertices[0];
        auto first_half = [&](Point<T> v) {
            T c = cross(anchor, v);
            return c > T(0) || (c == T(0) && dot(anchor, v) > T(0));
        };
        const bool tangent_first_half = first_half(tangent);
        Int lo = 0, hi = n;
        while (lo < hi) {
            Int mid = lo + (hi - lo) / 2;
            Point<T> edge = vertices[(mid + 1) % n] - vertices[mid];
            bool edge_first_half = first_half(edge);
            bool before_tangent = edge_first_half != tangent_first_half
                ? edge_first_half : cross(edge, tangent) > T(0);
            if (before_tangent)
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo == n ? 0 : lo;
    }
    // dot(direction, p) を最大にする頂点を返す。頂点のIDも保持する。条件・計算量・整数幅はsupport_indexと同じ。
    Point<T, Id> support_point(Point<T> direction) const {
        return (*this)[support_index(direction)];
    }
    // 直径の両端の、このConvexPolygon内の添字(i<=j)を返す。1点なら{0,0}、同距離なら任意の1組。
    // 添字は正規化後の頂点列に対応する。空の多角形は不可。回転calipersで時間 O(N)、追加領域 O(1)。
    // 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
    std::pair<Int, Int> diameter_pair() const {
        const Base &vertices = *this;
        Int n = this->size();
        assert(n > 0);
        T best = T(0);
        std::pair<Int, Int> indices{0, n == 1 ? 0 : 1};
        auto consider = [&](Int i, Int j) {
            T distance = norm(vertices[i] - vertices[j]);
            if (distance > best) {
                best = distance;
                indices = {std::min(i, j), std::max(i, j)};
            }
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
        return indices;
    }
    // 二乗直径を返す。空の多角形は不可。時間 O(N)、追加領域 O(1)。
    // 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
    T diameter() const {
        auto [i, j] = diameter_pair();
        return norm((*this)[i] - (*this)[j]);
    }
};
// 凸頂点列を前処理する。時間・領域 O(N)。
// 入力条件: vは凸多角形の周回順（時計回り・反時計回りどちらも可）の頂点列。共線・重複点も許す。
// 凸包構築・凸性検査は行わない。任意の点集合Pには initConvexPolygon(convex_hull(P)) を使う。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> auto initConvexPolygon(const std::vector<Point<T, Id>> &v) {
    return ConvexPolygon<T, Id>(v);
}
// 凸Polygonを前処理する。時間・領域 O(N)。
// 入力条件: polyは凸多角形の周回順（時計回り・反時計回りどちらも可）の頂点列。共線・重複点も許す。
// 凸包構築・凸性検査は行わない。任意の点集合Pには initConvexPolygon(convex_hull(P)) を使う。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> auto initConvexPolygon(const Polygon<T, Id> &poly) {
    return ConvexPolygon<T, Id>(poly);
}
// 正規化済みの凸多角形 a, b の頂点を合わせた凸包を返す。空・点・線分も可。
// 出力頂点が両方にある場合、IDはaのものを採用する。時間・領域 O(|a|+|b|)。
// 整数幅: B=8M^2（座標差の外積）。M=10^9ならT=long longで十分。
template <class T, class Id = void>
ConvexPolygon<T, Id> convex_merge(const ConvexPolygon<T, Id> &a, const ConvexPolygon<T, Id> &b) {
    if (a.empty())
        return b;
    if (b.empty())
        return a;

    using Vertex = Point<T, Id>;
    // 辞書順最小点から最大点まで、周回順・逆順のどちらでも辞書順に増加する。
    struct SortedVertices {
        const ConvexPolygon<T, Id> &vertices;
        std::size_t first = 0, last = vertices.size();

        bool empty() const { return first == last; }
        std::size_t next_index() const {
            return vertices[last - 1] < vertices[first] ? last - 1 : first;
        }
        const Vertex &front() const { return vertices[next_index()]; }
        const Vertex &pop() {
            const std::size_t index = next_index();
            if (index == first)
                ++first;
            else
                --last;
            return vertices[index];
        }
    };
    SortedVertices left{a}, right{b};
    std::vector<Vertex> points;
    points.reserve(a.size() + b.size());
    while (!left.empty() && !right.empty()) {
        auto &next = right.front() < left.front() ? right : left;
        points.push_back(next.pop());
    }
    while (!left.empty())
        points.push_back(left.pop());
    while (!right.empty())
        points.push_back(right.pop());

    ConvexPolygon<T, Id> result;
    result = detail::convexHullFromSorted(std::move(points), true);
    return result;
}
// 凸多角形の包含判定。時間 O(log(N+1))。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> bool contains(const ConvexPolygon<T, Id> &poly, Point<T> p, bool strict = false) {
    return poly.contains(p, strict);
}
// 凸多角形の境界判定。時間 O(log(N+1))。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> bool on_edge(const ConvexPolygon<T, Id> &poly, Point<T> p) {
    return poly.on_edge(p);
}
// 構築時に正規化した頂点列で、dot(direction, p) を最大にする頂点の添字を返す。
// 条件・計算量・整数幅はConvexPolygon::support_indexと同じ。
template <class T, class Id = void> Int support_index(const ConvexPolygon<T, Id> &poly, Point<T> direction) {
    return poly.support_index(direction);
}
// dot(direction, p) を最大にする頂点を返す。頂点のIDも保持する。条件・計算量・整数幅はsupport_indexと同じ。
template <class T, class Id = void> Point<T, Id> support_point(const ConvexPolygon<T, Id> &poly, Point<T> direction) {
    return poly.support_point(direction);
}
// 直径の両端のpoly内の添字(i<=j)を返す。1点なら{0,0}、同距離なら任意の1組。
// 空の多角形は不可。時間 O(N)、追加領域 O(1)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> std::pair<Int, Int> diameter_pair(const ConvexPolygon<T, Id> &poly) {
    return poly.diameter_pair();
}
// 二乗直径を返す。空の多角形は不可。時間 O(N)、追加領域 O(1)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> T diameter(const ConvexPolygon<T, Id> &poly) {
    return poly.diameter();
}
// 凸Polygonの二乗直径を返す。空の多角形は不可。時間・領域 O(N)。
// 入力条件: polyは凸多角形の周回順の頂点列。任意の点集合Pには diameter(convex_hull(P)) を使う。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T, class Id = void> T diameter(const Polygon<T, Id> &poly) {
    return initConvexPolygon(poly).diameter();
}
// 正規化した頂点列をコピーする。時間・領域 O(N)。
// 整数幅: 座標のコピーのみ。B=M、M=10^9ならint可。
template <class T, class Id = void> auto toPolygon(const ConvexPolygon<T, Id> &poly) {
    return poly.toPolygon();
}
namespace detail {
template <class T, class Id>
Polygon<T> convexCut(const std::vector<Point<T, Id>> &poly, Line<T> line) {
    assert(line.s != line.t);
    Polygon<T> result;
    auto append = [&](Point<T> p) {
        if (result.empty() || result.back() != p)
            result.push_back(p);
    };
    const Int n = static_cast<Int>(poly.size());
    for (Int i = 0; i < n; ++i) {
        const auto &a = poly[i], &b = poly[(i + 1) % n];
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
}
// 有向直線の左側と境界を座標だけの多角形として返す。浮動小数点またはFractionを使う。時間・出力領域 O(N)。
template <class T, class Id = void> requires (!std::is_integral_v<T>)
Polygon<T> convex_cut(const Polygon<T, Id> &poly, Line<T> line) {
    return detail::convexCut(poly, line);
}
// 前処理済み凸多角形を切断する。時間・領域 O(N)。
template <class T, class Id = void> requires (!std::is_integral_v<T>)
auto convex_cut(const ConvexPolygon<T, Id> &poly, Line<T> line) {
    return detail::convexCut(poly, line);
}
}
