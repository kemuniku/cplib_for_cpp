#pragma once
#include <cplib/geometry/convex_polygon.hpp>
#include <cplib/geometry/private/convex_arithmetic.hpp>
#include <optional>

namespace cplib {

template <class F> struct ConvexBoundingRectangle {
    std::array<Point<F>, 4> corners; // 反時計回り。点・線分では重複する。
    F area;
};

namespace detail {

template <class T> using ConvexMeasure = std::conditional_t<FloatingPoint<T>, T, long double>;

template <class F, class T> F convexMeasureValue(const T &value) {
    if constexpr (std::is_same_v<T, BigInt>)
        return std::stold(value.str());
    else
        return F(value);
}

template <bool Rectangle, class T, class Id>
auto convexCalipers(const ConvexPolygon<T, Id> &poly) {
    static_assert(FloatingPoint<T> || (convexInteger<T> && !std::is_same_v<T, BigInt>));
    using F = ConvexMeasure<T>;
    using Float = std::conditional_t<std::is_same_v<T, float> || std::is_same_v<T, double>, long double, T>;
    using W = std::conditional_t<FloatingPoint<T>, Float, BigInt>;
    using Approx = std::conditional_t<FloatingPoint<T>, Float, long double>;
    using Result = std::conditional_t<Rectangle, ConvexBoundingRectangle<F>, F>;
    const Int n = poly.size();
    if (n == 0)
        return std::optional<Result>{};
    auto approximate = [](const W &x) { return convexMeasureValue<Approx>(x); };
    auto approximate_point = [&](const Point<W> &p) -> Point<F> {
        return {F(approximate(p.x)), F(approximate(p.y))};
    };
    auto point = [&](Int i) { return convexPoint<W>(poly[i % n]); };
    auto origin = point(0);
    if (n <= 2) {
        if constexpr (Rectangle) {
            auto a = approximate_point(origin), b = approximate_point(point(n - 1));
            return std::optional<Result>(Result{{a, b, b, a}, F(0)});
        } else {
            return std::optional<Result>(F(0));
        }
    }

    auto vertex = [&](Int i) { return point(i) - origin; };
    auto first_dir = vertex(1);
    Int top = convexSupportIndex(n, vertex, convexNormal(first_dir));
    Int right = 0, left = 0;
    if constexpr (Rectangle) {
        right = convexSupportIndex(n, vertex, first_dir);
        left = convexSupportIndex(n, vertex, -first_dir);
    }
    auto advance = [&](Int &index, const Point<W> &d) {
        while (dot(d, vertex(index + 1) - vertex(index)) > W(0))
            index = (index + 1) % n;
    };
    W best_numerator(0), best_denominator(0), best_left(0), best_right(0), best_top(0), best_bottom(0);
    Point<W> best_direction;
    for (Int i = 0; i < n; ++i) {
        auto p = vertex(i);
        auto d = vertex(i + 1) - p;
        auto normal = convexNormal(d);
        advance(top, normal);
        W bottom = dot(normal, p), top_projection = dot(normal, vertex(top));
        W height = top_projection - bottom;
        W d2 = norm(d), numerator = height * height;
        W l(0), r(0);
        if constexpr (Rectangle) {
            advance(right, d);
            advance(left, -d);
            l = dot(d, vertex(left));
            r = dot(d, vertex(right));
            numerator = height * (r - l);
        }
        if (i != 0 && !(numerator * best_denominator < best_numerator * d2))
            continue;
        best_numerator = std::move(numerator);
        best_denominator = std::move(d2);
        if constexpr (Rectangle) {
            best_direction = std::move(d);
            best_left = std::move(l);
            best_right = std::move(r);
            best_bottom = std::move(bottom);
            best_top = std::move(top_projection);
        }
    }

    Approx denominator = approximate(best_denominator);
    Approx measure = approximate(best_numerator) / denominator;
    if constexpr (Rectangle) {
        auto normal = convexNormal(best_direction);
        auto corner = [&](const W &x, const W &y) {
            auto offset = best_direction * x + normal * y;
            return Point<F>{F(approximate(origin.x) + approximate(offset.x) / denominator),
                            F(approximate(origin.y) + approximate(offset.y) / denominator)};
        };
        return std::optional<Result>(Result{{corner(best_left, best_bottom), corner(best_right, best_bottom),
                                             corner(best_right, best_top), corner(best_left, best_top)}, F(measure)});
    } else {
        return std::optional<Result>(F(cplib::sqrt(measure)));
    }
}

}

// 正規化済み凸多角形の平行支持線間の最小距離。空は nullopt、点・線分は 0。
// 回転キャリパー法で O(N) 回の算術・比較、追加領域 O(1)。整数入力は long double で返す。
template <class T, class Id>
std::optional<detail::ConvexMeasure<T>> convex_minimum_width(const ConvexPolygon<T, Id> &poly) {
    return detail::convexCalipers<false>(poly);
}

// 正規化済み凸多角形を囲む最小面積の長方形。空は nullopt。点・線分では面積 0。
// 回転キャリパー法で O(N) 回の算術・比較、追加領域 O(1)。整数入力は長方形を long double で返す。
template <class T, class Id>
std::optional<ConvexBoundingRectangle<detail::ConvexMeasure<T>>>
convex_minimum_area_rectangle(const ConvexPolygon<T, Id> &poly) {
    return detail::convexCalipers<true>(poly);
}

}
