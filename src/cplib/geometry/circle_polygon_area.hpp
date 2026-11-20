#pragma once
#include <cplib/geometry/circle.hpp>
#include <cplib/geometry/polygon.hpp>

namespace cplib::detail {
// 原点と有向辺a->bの三角形と、原点中心の円盤との符号付き共通面積。
template <class F> F circlePolygonEdgeArea(Point<F> a, Point<F> b, F radius) {
    auto triangle = [](Point<F> p, Point<F> q) {
        F scale = std::max({abs(p.x), abs(p.y), abs(q.x), abs(q.y)});
        if (scale == 0)
            return F(0);
        return (cplib::cross(p / scale, q / scale) * scale / 2) * scale;
    };
    auto sector = [&](Point<F> p, Point<F> q) {
        F lp = circleLength(p), lq = circleLength(q);
        if (lp == 0 || lq == 0)
            return F(0);
        p /= lp;
        q /= lq;
        return (atan2(cplib::cross(p, q), cplib::dot(p, q)) * radius / 2) * radius;
    };
    if (exact_equal(a, b))
        return 0;
    if (circleLength(a) <= radius && circleLength(b) <= radius)
        return triangle(a, b);

    auto direction = b - a;
    auto unit = direction / circleLength(direction);
    // 交点をa+t*(b-a)で作ると、遠い端点からの補間で桁落ちする。
    // 垂線の足と、辺方向の距離で円周上の交点を直接構築する。
    F height = cplib::cross(unit, a), distance = abs(height);
    if (distance >= radius)
        return sector(a, b);
    F ratio = distance / radius;
    F width = radius * sqrt((1 - ratio) * (1 + ratio));
    F start = cplib::dot(a, unit), finish = cplib::dot(b, unit);
    if (finish <= -width || start >= width)
        return sector(a, b);
    auto foot = circlePerpendicular(unit) * height;
    auto entry = start >= -width ? a : foot - unit * width;
    auto exit = finish <= width ? b : foot + unit * width;
    return sector(a, entry) + triangle(entry, exit) + sector(exit, b);
}
}

namespace cplib {
// 円盤と単純多角形の共通部分面積を非負の近似値で返す。凹多角形・両方の頂点順に対応。
// 3頂点未満または半径0なら0。重複頂点・始点の再掲も許す。自己交差は対象外。
// 戻り値はCircleResultScalar<T,S>、内部計算はlong double以上の精度。整数・Fractionも近似する。
// long doubleより高精度な入力型は、内部計算でもその精度を保つ。
// 座標・差分・面積は使用する浮動小数点型で有限に表現できること。
// 時間 O(N) 算術操作と近似変換費用、追加領域 O(1)。厳密円の中間値は指定基底型に収まること。
template <class T, class S, class Id = void>
// 整数幅: 厳密円の変換にcircle.hppのh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
// 多角形の座標は積・差を取る前に浮動小数点へ変換するので、整数の追加幅は不要。
CircleResultScalar<T, S> intersection_area(const Circle<T> &circle, const Polygon<S, Id> &polygon) {
    using F = CircleResultScalar<T, S>;
    using Work = std::common_type_t<long double, F>;
    if (polygon.size() < 3)
        return F(0);
    Point<Work> center;
    Work radius;
    if constexpr (cplib::isFloatingPoint<T>) {
        center = detail::circleFloatPointAs<Work>(circle.center());
        radius = circle.radius();
    } else {
        auto data = detail::circleExactData(circle);
        center = {detail::circleApprox<Work>(data.first.x),
                  detail::circleApprox<Work>(data.first.y)};
        radius = detail::circleApprox<Work>(data.second, true);
    }
    if (radius == 0)
        return F(0);
    auto relative = [&](const Point<S> &point) {
        auto p = detail::circleFloatPointAs<Work>(point) - center;
        detail::circleCheckPoint(p);
        return p;
    };
    Work sum = 0, correction = 0;
    auto previous = relative(polygon.back());
    for (const auto &point : polygon) {
        auto current = relative(point);
        Work term = detail::circlePolygonEdgeArea(previous, current, radius);
        // 円の外側の辺による扇形の相殺などで、加算誤差が蓄積するのを抑える。
        Work next = sum + term;
        correction += abs(sum) >= abs(term)
            ? (sum - next) + term : (term - next) + sum;
        sum = next;
        previous = current;
    }
    F result = static_cast<F>(abs(sum + correction));
    if (!isfinite(result))
        throw std::invalid_argument("円と多角形の共通部分面積は有限で表現できる必要があります");
    return result;
}

// intersection_area(circle, polygon)と同じ。時間 O(N) 算術操作と近似変換費用、追加領域 O(1)。
template <class S, class T, class Id = void>
// 整数幅: 引数順を逆にした版と同じh+1ビット。M=10^9の一般の整数3点円ではBigInt等。
CircleResultScalar<T, S> intersection_area(const Polygon<S, Id> &polygon, const Circle<T> &circle) {
    return cplib::intersection_area(circle, polygon);
}
}
