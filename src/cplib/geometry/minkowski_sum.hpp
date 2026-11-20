#pragma once
#include <cplib/geometry/convex_polygon.hpp>
#include <cplib/geometry/private/convex_arithmetic.hpp>
#include <stdexcept>
#include <tuple>

namespace cplib {

// 正規化済み凸多角形のミンコフスキー和。空・点・線分にも対応する。
// 結果の頂点は座標のみ。O(N+M) 回の算術・比較、追加領域 O(1)（出力を除く）。
// 固定幅整数では和の座標が T に収まらない場合に overflow_error。
template <class T, class A, class B>
ConvexPolygon<T> minkowski_sum(const ConvexPolygon<T, A> &a, const ConvexPolygon<T, B> &b) {
    using W = detail::ConvexWide<T>;
    ConvexPolygon<T> result;
    if (a.empty() || b.empty())
        return result;

    auto start = [](const auto &p) {
        return std::min_element(p.begin(), p.end(), [](const auto &l, const auto &r) {
            return std::tie(l.y, l.x) < std::tie(r.y, r.x);
        }) - p.begin();
    };
    auto edge = [](const auto &p, std::size_t i) {
        return detail::convexPoint<W>(p[(i + 1) % p.size()]) - detail::convexPoint<W>(p[i]);
    };
    auto half = [](const Point<W> &v) {
        return v.y < W(0) || (v.y == W(0) && v.x < W(0));
    };
    auto compare = [&](const Point<W> &u, const Point<W> &v) {
        if (half(u) != half(v))
            return half(u) ? 1 : -1;
        auto turn = cross(u, v);
        return (turn < W(0)) - (turn > W(0));
    };
    auto add = [](const T &x, const T &y) {
        if constexpr (detail::convexInteger<T> && !std::is_same_v<T, BigInt>) {
            T sum;
            if (__builtin_add_overflow(x, y, &sum))
                throw std::overflow_error("minkowski_sum coordinate is out of range");
            return sum;
        } else {
            return x + y;
        }
    };

    std::size_t i = start(a), j = start(b);
    std::size_t remaining_a = a.size() > 1 ? a.size() : 0;
    std::size_t remaining_b = b.size() > 1 ? b.size() : 0;
    result.reserve(std::max<std::size_t>(1, remaining_a + remaining_b));
    do {
        result.push_back({add(a[i].x, b[j].x), add(a[i].y, b[j].y)});
        int order = remaining_a == 0 ? 1 : remaining_b == 0 ? -1 : compare(edge(a, i), edge(b, j));
        if (remaining_a != 0 && order <= 0) {
            --remaining_a;
            i = (i + 1) % a.size();
        }
        if (remaining_b != 0 && order >= 0) {
            --remaining_b;
            j = (j + 1) % b.size();
        }
    } while (remaining_a != 0 || remaining_b != 0);
    std::rotate(result.begin(), std::min_element(result.begin(), result.end()), result.end());
    return result;
}

}
