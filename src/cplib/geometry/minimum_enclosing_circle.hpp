#pragma once
#include <cplib/geometry/circle.hpp>
#include <random>

namespace cplib::detail {
template <class F> Circle<F> mecFloatDiameter(Point<F> a, Point<F> b) {
    auto center = a + (b - a) * 0.5;
    return Circle<F>(center, std::max(circleLength(a - center), circleLength(b - center)));
}
template <class F> Circle<F> mecFloatThree(Point<F> a, Point<F> b, Point<F> c) {
    auto u = b - a, v = c - a;
    F scale = std::max({abs(u.x), abs(u.y), abs(v.x), abs(v.y)});
    if (scale == 0)
        return Circle<F>(a, 0);
    auto un = u / scale, vn = v / scale;
    F determinant = 2 * cplib::cross(un, vn);
    if (determinant == 0) {
        auto result = mecFloatDiameter(a, b);
        for (auto candidate : {mecFloatDiameter(a, c), mecFloatDiameter(b, c)})
            if (candidate.radius() > result.radius())
                result = candidate;
        return result;
    }
    Point<F> offset{
        (vn.y * cplib::norm(un) - un.y * cplib::norm(vn)) / determinant * scale,
        (un.x * cplib::norm(vn) - vn.x * cplib::norm(un)) / determinant * scale};
    auto center = a + offset;
    return Circle<F>(center, std::max({circleLength(a - center), circleLength(b - center),
                                           circleLength(c - center)}));
}
template <class T> bool mecRadiusLess(const Circle<T> &a, const Circle<T> &b) {
    auto x = a.radius_squared_exact(), y = b.radius_squared_exact();
    return x.num * y.den < y.num * x.den;
}
template <class T> auto mecExactThree(const Point<T> &a, const Point<T> &b, const Point<T> &c) {
    try {
        return initCircle(a, b, c);
    } catch (const std::invalid_argument &) {
        auto result = initDiameterCircle(a, b);
        for (auto candidate : {initDiameterCircle(a, c), initDiameterCircle(b, c)})
            if (mecRadiusLess(result, candidate))
                result = std::move(candidate);
        return result;
    }
}
}

namespace cplib {
// 点集合の最小包含円を返す。期待 O(N)、最悪 O(N³) 算術操作、追加領域 O(N)。
// 非floatは全判定を指定基底型で厳密に行い、中間値がその型に収まることを前提とする。
// floatは平行移動・尺度正規化・最後の半径補正を行う。多倍長演算の費用は桁数に依存する。
// 空入力はinvalid_argument。入力とグローバル乱数状態を変えず、同じ入力とseedで再現できる。
template <class T, class Id = void>
// 整数幅: 入力座標上限Mから保持座標は4M以下、尺度は1または2、classifyのHも4M以下。
// k=ceil(log2(4M+1)), h=6k+11として、max(W(96*(4M)^4),2h+1)ビットで十分。
// 2h+1は共線時に候補円の半径二乗を分子・分母の交差積で比較するため。Nによる幅の増加はない。
// M=10^9ならh=203、407ビットで十分。__int128では保証できず、入力からBigInt等を使う。
// Fraction入力はcircle.hppの規則で保持点構築とclassifyの幅も別途見積もる。
auto minimum_enclosing_circle(const std::vector<Point<T, Id>> &points, Int seed = 0, long double tolerance = 1e-12L) {
    if (points.empty())
        throw std::invalid_argument("最小包含円の入力は非空である必要があります");
    std::mt19937_64 rng(static_cast<UInt>(seed));
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = detail::CircleFloatType<T>;
        detail::circleCheckTolerance(tolerance);
        std::vector<Point<F>> original;
        original.reserve(points.size());
        for (const auto &point : points) {
            auto p = detail::circleFloatPointAs<F>(point);
            detail::circleCheckPoint(p);
            original.push_back(p);
        }
        auto origin = original[0];
        F scale = 0;
        for (auto point : original)
            scale = std::max(scale, detail::circleLength(point - origin));
        if (scale == 0)
            return Circle<F>(origin, 0);
        std::vector<Point<F>> shuffled;
        shuffled.reserve(points.size());
        for (auto point : original)
            shuffled.push_back((point - origin) / scale);
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        Circle<F> circle(shuffled[0], 0);
        for (std::size_t i = 0; i < shuffled.size(); ++i) {
            auto p = shuffled[i];
            if (circle.contains(p, tolerance))
                continue;
            circle = Circle<F>(p, 0);
            for (std::size_t j = 0; j < i; ++j) {
                auto q = shuffled[j];
                if (circle.contains(q, tolerance))
                    continue;
                circle = detail::mecFloatDiameter(p, q);
                for (std::size_t k = 0; k < j; ++k)
                    if (!circle.contains(shuffled[k], tolerance))
                        circle = detail::mecFloatThree(p, q, shuffled[k]);
            }
        }
        auto center = origin + circle.center() * scale;
        F radius = 0;
        for (auto point : original)
            radius = std::max(radius, detail::circleLength(point - center));
        return Circle<F>(center, radius);
    } else {
        std::vector<Point<T>> shuffled(points.begin(), points.end());
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        auto circle = initCircle(shuffled[0], shuffled[0], shuffled[0]);
        for (std::size_t i = 0; i < shuffled.size(); ++i) {
            const auto &p = shuffled[i];
            if (circle.contains(p))
                continue;
            circle = initCircle(p, p, p);
            for (std::size_t j = 0; j < i; ++j) {
                const auto &q = shuffled[j];
                if (circle.contains(q))
                    continue;
                circle = initDiameterCircle(p, q);
                for (std::size_t k = 0; k < j; ++k)
                    if (!circle.contains(shuffled[k]))
                        circle = detail::mecExactThree(p, q, shuffled[k]);
            }
        }
        return circle;
    }
}
}
