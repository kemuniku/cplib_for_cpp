#pragma once
#include <cplib/geometry/base.hpp>

namespace cplib {
inline constexpr Int ANGLE_0 = 0, ANGLE_0_90 = 1, ANGLE_90 = 2, ANGLE_90_180 = 3, ANGLE_180 = 4,
                     ANGLE_180_270 = -3, ANGLE_270 = -2, ANGLE_270_360 = -1;

// 非零ベクトル p, q のなす角を八方位で返す。計算量 O(1)。
// 整数幅: ベクトル成分をM以下としてB=2M^2。M=10^9なら62ビット、long long。
template <class T> Int angle(Point<T> p, Point<T> q) {
    assert(!(geometry_eq(p.x, 0) && geometry_eq(p.y, 0)) &&
           !(geometry_eq(q.x, 0) && geometry_eq(q.y, 0)));
    T d = dot(p, q), c = cross(p, q);
    if (geometry_eq(c, 0))
        return geometry_gt(d, 0) ? ANGLE_0 : ANGLE_180;
    if (geometry_eq(d, 0))
        return geometry_gt(c, 0) ? ANGLE_90 : ANGLE_270;
    if (geometry_gt(d, 0) && geometry_gt(c, 0))
        return ANGLE_0_90;
    if (geometry_lt(d, 0) && geometry_gt(c, 0))
        return ANGLE_90_180;
    if (geometry_lt(d, 0) && geometry_lt(c, 0))
        return ANGLE_180_270;
    if (geometry_gt(d, 0) && geometry_lt(c, 0))
        return ANGLE_270_360;
    return 0;
}

// 2直線の方向ベクトルのなす角を八方位で返す。計算量 O(1)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T> Int angle(const Line<T> &a, const Line<T> &b) {
    return angle(toVector(a), toVector(b));
}

// 2線分の方向ベクトルのなす角を八方位で返す。計算量 O(1)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T> Int angle(const Segment<T> &a, const Segment<T> &b) {
    return angle(toVector(a), toVector(b));
}

// 直線と線分の方向ベクトルのなす角を八方位で返す。計算量 O(1)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T> Int angle(const Line<T> &a, const Segment<T> &b) {
    return angle(a, toLine(b));
}

// 線分と直線の方向ベクトルのなす角を八方位で返す。計算量 O(1)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T> Int angle(const Segment<T> &a, const Line<T> &b) {
    return angle(toLine(a), b);
}

// a, bが平行かどうかを判定する。計算量 O(1)。
// 整数幅: angleと同じ。Point同士はB=2M^2、直線・線分同士はB=8M^2。M=10^9ならlong long。
template <class P, class Q> bool is_parallel(const P &a, const Q &b) {
    Int d = angle(a, b);
    return d == ANGLE_0 || d == ANGLE_180;
}

// a, bが直角かどうかを判定する。計算量 O(1)。
// 整数幅: angleと同じ。Point同士はB=2M^2、直線・線分同士はB=8M^2。M=10^9ならlong long。
template <class P, class Q> bool is_orthogonal(const P &a, const Q &b) {
    Int d = angle(a, b);
    return d == ANGLE_90 || d == ANGLE_270;
}
}
