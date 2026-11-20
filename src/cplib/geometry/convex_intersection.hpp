#pragma once
#include <cplib/geometry/convex_polygon_queries.hpp>
#include <cplib/geometry/half_plane_intersection.hpp>

namespace cplib {

// 正規化済み凸多角形の共通部分。空・点・線分も返す。
// 整数座標の結果は厳密な Fraction<BigInt>。浮動小数点では半平面交差と同じ精度・入力条件。
// O((N+M) log(N+M+1)) 回の体演算・比較、領域 O(N+M)。
template <class T, class A, class B>
ConvexPolygon<detail::ConvexField<T>> convex_intersection(const ConvexPolygon<T, A> &a,
                                                        const ConvexPolygon<T, B> &b) {
    using R = detail::ConvexField<T>;
    using P = Point<R>;
    if (a.empty() || b.empty())
        return {};

    auto clip_small = [](const auto &small, const auto &poly) -> ConvexPolygon<R> {
        if (small.size() == 1) {
            if (convex_tangents(poly, Point<T>(small.front())))
                return {};
            return {detail::convexPoint<R>(small.front())};
        }

        auto section = convex_line_intersection(poly, Line<T>{small.front(), small.back()});
        if (!section)
            return {};
        // 正規化済みの線分は辞書順。section の両端も同じ順に並ぶ。
        P first = std::max(section->s, detail::convexPoint<R>(small.front()));
        P last = std::min(section->t, detail::convexPoint<R>(small.back()));
        if (last < first)
            return {};
        if (exact_equal(first, last))
            return {first};
        return {first, last};
    };
    if (a.size() <= 2)
        return clip_small(a, b);
    if (b.size() <= 2)
        return clip_small(b, a);

    auto for_edges = [&](const auto &emit) {
        auto append = [&](const auto &poly) {
            for (std::size_t i = 0; i < poly.size(); ++i)
                emit(poly[i], poly[(i + 1) % poly.size()]);
        };
        append(a);
        append(b);
    };
    if constexpr (detail::convexInteger<T>) {
        std::vector<std::array<BigInt, 3>> constraints;
        constraints.reserve(a.size() + b.size());
        for_edges([&](const auto &u, const auto &v) {
            auto p = detail::convexPoint<BigInt>(u), q = detail::convexPoint<BigInt>(v);
            auto normal = detail::convexNormal(q - p);
            constraints.push_back({-normal.x, -normal.y, dot(normal, p)});
        });
        return half_plane_intersection(constraints).toConvexPolygon();
    } else {
        std::vector<Line<R>> lines;
        lines.reserve(a.size() + b.size());
        for_edges([&](const auto &p, const auto &q) { lines.push_back({p, q}); });
        return half_plane_intersection(lines).toConvexPolygon();
    }
}

// 共通部分の面積。整数座標の入力でも面積を約分可能な有理数で返す。
template <class T, class A, class B>
detail::ConvexField<T> convex_intersection_area(const ConvexPolygon<T, A> &a, const ConvexPolygon<T, B> &b) {
    auto intersection = convex_intersection(a, b);
    if (intersection.size() < 3)
        return detail::ConvexField<T>(0);
    return area(intersection.toPolygon());
}

}
