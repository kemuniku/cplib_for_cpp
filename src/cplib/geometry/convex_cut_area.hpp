#pragma once
#include <cplib/geometry/convex_polygon.hpp>

namespace cplib {
// 正規化済みの凸多角形をコピーし、有向直線の左側に残る面積を高速に求める。
// 構築は時間・領域 O(N)、各問い合わせは O(log N)。
template <class T> requires (!std::is_integral_v<T>)
class ConvexCutArea {
    std::vector<Point<T>> vertices_;
    std::vector<T> prefix_;

    // 辺の方向を、最初の辺を基準に偏角順で比較する。
    Int support_index(Point<T> direction) const {
        const Int n = static_cast<Int>(vertices_.size());
        Point<T> anchor = vertices_[1] - vertices_[0];
        auto half = [&](Point<T> v) {
            T c = cross(anchor, v);
            return c > T(0) || (c == T(0) && dot(anchor, v) > T(0)) ? 0 : 1;
        };
        auto less_angle = [&](Point<T> a, Point<T> b) {
            int ha = half(a), hb = half(b);
            return ha != hb ? ha < hb : cross(a, b) > T(0);
        };
        Int lo = 0, hi = n;
        while (lo < hi) {
            Int mid = lo + (hi - lo) / 2;
            Point<T> edge = vertices_[(mid + 1) % n] - vertices_[mid];
            if (less_angle(edge, direction))
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo == n ? 0 : lo;
    }

public:
    explicit ConvexCutArea(const ConvexPolygon<T> &poly) : vertices_(poly.begin(), poly.end()) {
        const Int n = static_cast<Int>(vertices_.size());
        prefix_.reserve(2 * n + 1);
        prefix_.push_back(T(0));
        for (Int i = 0; i < 2 * n; ++i)
            prefix_.push_back(prefix_.back() + cross(vertices_[i % n], vertices_[(i + 1) % n]));
    }

    // line の左側と境界に残る部分の面積を返す。計算量 O(log N)。
    T area_left(const Line<T> &line) const {
        assert(line.s != line.t);
        const Int n = static_cast<Int>(vertices_.size());
        if (n < 3)
            return T(0);
        Point<T> direction = toVector(line);
        auto side = [&](Int i) { return cross(direction, vertices_[i % n] - line.s); };
        Int minimum = support_index(direction);
        Int maximum = support_index(-direction);
        Int min_sign = detail::convexSign(side(minimum));
        Int max_sign = detail::convexSign(side(maximum));
        if (min_sign >= 0)
            return prefix_[n] / T(2);
        if (max_sign <= 0)
            return T(0);

        if (maximum <= minimum)
            maximum += n;
        // 最小点から最大点へ進む辺上で、初めて左側に入る頂点を探す。
        Int lo = minimum + 1, hi = maximum;
        while (lo < hi) {
            Int mid = lo + (hi - lo) / 2;
            if (detail::convexSign(side(mid)) >= 0)
                hi = mid;
            else
                lo = mid + 1;
        }
        Int first = lo;
        // 最大点から最小点へ進む辺上で、初めて左側を出る頂点を探す。
        lo = maximum + 1;
        hi = minimum + n;
        while (lo < hi) {
            Int mid = lo + (hi - lo) / 2;
            if (detail::convexSign(side(mid)) <= 0)
                hi = mid;
            else
                lo = mid + 1;
        }
        Int after_last = lo;
        auto intersection = [&](Int a, Int b) {
            Point<T> p = vertices_[a % n], q = vertices_[b % n];
            T fp = side(a), fq = side(b);
            if (detail::convexSign(fp) == 0)
                return p;
            if (detail::convexSign(fq) == 0)
                return q;
            return p + (q - p) * (fp / (fp - fq));
        };
        Point<T> entry = intersection(first - 1, first);
        Point<T> exit = intersection(after_last - 1, after_last);
        Int last = after_last - 1;
        T twice_area = cross(entry, vertices_[first % n]) +
                       (prefix_[last] - prefix_[first]) +
                       cross(vertices_[last % n], exit) + cross(exit, entry);
        T result = twice_area / T(2);
        return result > T(0) ? result : T(0);
    }
};

template <class T> requires (!std::is_integral_v<T>)
ConvexCutArea<T> initConvexCutArea(const ConvexPolygon<T> &poly) {
    return ConvexCutArea<T>(poly);
}
}
