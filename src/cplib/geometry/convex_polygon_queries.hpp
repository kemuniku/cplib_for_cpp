#pragma once
#include <cplib/geometry/convex_polygon.hpp>
#include <cplib/geometry/private/convex_arithmetic.hpp>
#include <optional>

namespace cplib {

// 多角形全体が、外部点 p から接点への有向直線の左側、右側にある頂点の添字を順に返す。
// 空、内部・境界上の p では nullopt。接触する辺がある場合はその端点のいずれかを返す。
// 正規化済み凸多角形について O(log(N+1)) 回の算術・比較、追加領域 O(1)。
template <class T, class Id>
std::optional<std::pair<Int, Int>> convex_tangents(const ConvexPolygon<T, Id> &poly, Point<T> p) {
    using W = detail::ConvexWide<T>;
    const Int n = poly.size();
    if (n == 0)
        return std::nullopt;
    auto vertex = [&](Int i) { return detail::convexPoint<W>(poly[i % n]); };
    const Point<W> q = detail::convexPoint<W>(p);
    auto side = [&](Int i) { return cross(vertex(i + 1) - vertex(i), q - vertex(i)); };
    const W zero(0);
    if (n == 1)
        return exact_equal(vertex(0), q) ? std::nullopt : std::optional(std::pair<Int, Int>{0, 0});
    if (n == 2) {
        auto turn = side(0);
        if (turn == zero && dot(q - vertex(0), q - vertex(1)) <= zero)
            return std::nullopt;
        return turn < zero ? std::pair<Int, Int>{1, 0} : std::pair<Int, Int>{0, 1};
    }

    Int visible = 0;
    if (side(0) >= zero) {
        if (side(n - 1) < zero) {
            visible = n - 1;
        } else {
            Int lo = 1, hi = n - 1;
            const auto origin = vertex(0);
            while (hi - lo > 1) {
                Int mid = (lo + hi) / 2;
                if (cross(vertex(mid) - origin, q - origin) >= zero)
                    lo = mid;
                else
                    hi = mid;
            }
            visible = lo;
            if (side(visible) >= zero) {
                // 扇の境界の延長上では、端点に隣接する辺も調べる。
                if (cross(vertex(1) - origin, q - origin) == zero && side(1) < zero)
                    visible = 1;
                else if (cross(vertex(n - 1) - origin, q - origin) == zero && side(n - 2) < zero)
                    visible = n - 2;
                else
                    return std::nullopt;
            }
        }
    }

    // 可視辺の内向き法線で最も遠い頂点を取る。そこに接する少なくとも一辺は非可視。
    auto normal = detail::convexNormal(vertex(visible + 1) - vertex(visible));
    Int far = detail::convexSupportIndex(n, vertex, normal);
    Int invisible = side(far) >= zero ? far : (far + n - 1) % n;
    if (invisible < visible)
        invisible += n;
    Int lo = visible, hi = invisible;
    while (hi - lo > 1) {
        Int mid = (lo + hi) / 2;
        if (side(mid) < zero)
            lo = mid;
        else
            hi = mid;
    }
    const Int left = (lo + 1) % n;
    lo = invisible;
    hi = visible + n;
    while (hi - lo > 1) {
        Int mid = (lo + hi) / 2;
        if (side(mid) >= zero)
            lo = mid;
        else
            hi = mid;
    }
    return std::pair<Int, Int>{left, hi % n};
}

// 正規化済み凸多角形と直線の共通部分。空は nullopt、接点は両端が同じ線分。
// 端点は直線の向きに沿って並べる。整数座標の結果は Fraction<BigInt>。
// O(log(N+1)) 回の算術・比較、追加領域 O(1)。直線の端点は相異なること。
template <class T, class Id>
auto convex_line_intersection(const ConvexPolygon<T, Id> &poly, const Line<T> &line)
    -> std::optional<Segment<detail::ConvexField<T>>> {
    using W = detail::ConvexWide<T>;
    using R = detail::ConvexField<T>;
    using P = Point<R>;
    const Int n = poly.size();
    if (n == 0)
        return std::nullopt;
    const auto start = detail::convexPoint<W>(line.s), finish = detail::convexPoint<W>(line.t);
    assert(!exact_equal(start, finish));
    auto vertex = [&](Int i) { return detail::convexPoint<W>(poly[i % n]); };
    auto lifted = [&](Int i) { return detail::convexPoint<R>(poly[i % n]); };
    auto direction = finish - start;
    auto normal = detail::convexNormal(direction);
    auto side = [&](Int i) { return dot(normal, vertex(i) - start); };
    auto ordered = [&](P a, P b) {
        auto d = detail::convexPoint<R>(line.t) - detail::convexPoint<R>(line.s);
        if (dot(d, b - a) < R(0))
            std::swap(a, b);
        return Segment<R>{std::move(a), std::move(b)};
    };
    Int low = detail::convexSupportIndex(n, vertex, -normal);
    Int high = detail::convexSupportIndex(n, vertex, normal);
    W lower = side(low), upper = side(high), zero(0);
    if (lower > zero || upper < zero)
        return std::nullopt;

    if (lower == zero || upper == zero) {
        if (lower == upper)
            return ordered(lifted(0), lifted(n - 1));
        Int index = lower == zero ? low : high;
        Int neighbor = index;
        if (side(index + 1) == zero)
            neighbor = index + 1;
        else if (side(index + n - 1) == zero)
            neighbor = index + n - 1;
        return ordered(lifted(index), lifted(neighbor));
    }

    auto intersection = [&](Int lo, Int hi, bool increasing) -> P {
        while (hi - lo > 1) {
            Int mid = (lo + hi) / 2;
            auto value = side(mid);
            if (value == zero)
                return lifted(mid);
            if ((value < zero) == increasing)
                lo = mid;
            else
                hi = mid;
        }
        auto a = lifted(lo), b = lifted(hi);
        auto rstart = detail::convexPoint<R>(line.s);
        auto rdir = detail::convexPoint<R>(line.t) - rstart;
        R fa = cross(rdir, a - rstart), fb = cross(rdir, b - rstart);
        return a + (b - a) * (fa / (fa - fb));
    };
    if (high < low)
        high += n;
    return ordered(intersection(low, high, true), intersection(high, low + n, false));
}

}
