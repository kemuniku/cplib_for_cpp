#pragma once
#include <cplib/geometry/base.hpp>

namespace cplib {
// 直線 l への点 p の射影。計算量 O(1)。整数座標では途中の除算が整数除算になる。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
// 除算後に方向ベクトルを掛ける実装なので、整数版の中間値は三乗にはならない。
template <class T> Point<T> projection(const Line<T> &l, Point<T> p) {
    T t = dot(p - l.s, toVector(l)) / norm(toVector(l));
    return l.s + toVector(l) * t;
}

// 点 p の直線 l に対する反射。計算量 O(1)。整数座標では射影の途中が整数除算になる。
// 整数幅: B=max(8M^2,12M)。射影と反射の加減算・2倍を含む。M=10^9なら64ビット、long long。
template <class T> Point<T> reflection(const Line<T> &l, Point<T> p) {
    return p + (projection(l, p) - p) * 2;
}

// 線分を延長した直線への点pの射影を返す。計算量 O(1)。整数座標では途中が整数除算になる。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
// 除算後に方向ベクトルを掛ける実装なので、整数版の中間値は三乗にはならない。
template <class T> auto projection(const Segment<T> &l, Point<T> p) {
    return projection(toLine(l), p);
}

// 線分を延長した直線に対する点pの反射を返す。計算量 O(1)。整数座標では途中が整数除算になる。
// 整数幅: B=max(8M^2,12M)。射影と反射の加減算・2倍を含む。M=10^9なら64ビット、long long。
template <class T> auto reflection(const Segment<T> &l, Point<T> p) {
    return reflection(toLine(l), p);
}
}
