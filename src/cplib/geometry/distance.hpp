#pragma once
#include <cplib/geometry/intersect.hpp>

namespace cplib {
// 点 a, b のユークリッド距離の2乗。計算量 O(1)。
// 整数幅: B=8M^2（座標差の内積・外積・二乗）。M=10^9なら64ビット、T=long longで十分。
template <class T> T norm(Point<T> a, Point<T> b) {
    return norm(a - b);
}

// 点 p と直線 l のユークリッド距離の2乗。計算量 O(1)。整数座標では整数除算で切り捨てる。
// 整数幅: B=64M^4。外積をさらに二乗してから割る。M=10^9なら127ビット、T=__int128が必要。
// 最終的な距離二乗がlong longに収まっても、中間のarea*areaは収まらない場合がある。
template <class T> T norm(Point<T> p, const Line<T> &l) {
    T area = cross(toVector(l), p - l.s);
    return area * area / norm(toVector(l));
}

// 点 p と線分 s のユークリッド距離の2乗。計算量 O(1)。整数座標では直線への距離に整数除算を用いる。
// 整数幅: B=64M^4。外積をさらに二乗してから割る。M=10^9なら127ビット、T=__int128が必要。
// 最終的な距離二乗がlong longに収まっても、中間のarea*areaは収まらない場合がある。
template <class T> T norm(Point<T> p, const Segment<T> &s) {
    if (geometry_lt(dot(toVector(s), p - s.s), 0))
        return norm(p - s.s);
    if (geometry_lt(dot(-toVector(s), p - s.t), 0))
        return norm(p - s.t);
    return norm(p, toLine(s));
}

// 線分 a, b のユークリッド距離の2乗。計算量 O(1)。整数座標では途中に整数除算を用いる。
// 整数幅: B=64M^4。外積をさらに二乗してから割る。M=10^9なら127ビット、T=__int128が必要。
// 最終的な距離二乗がlong longに収まっても、中間のarea*areaは収まらない場合がある。
template <class T> T norm(const Segment<T> &a, const Segment<T> &b) {
    if (intersect(a, b))
        return T(0);
    return std::min({norm(a.s, b), norm(a.t, b), norm(b.s, a), norm(b.t, a)});
}

// 点 a, b のユークリッド距離。計算量 O(1)。
template <cplib::FloatingPoint T> T distance(Point<T> a, Point<T> b) {
    return sqrt(norm(a, b));
}

// 点 p と直線 l のユークリッド距離。計算量 O(1)。
template <cplib::FloatingPoint T> T distance(Point<T> p, const Line<T> &l) {
    return sqrt(norm(p, l));
}

// 点 p と線分 s のユークリッド距離。計算量 O(1)。
template <cplib::FloatingPoint T> T distance(Point<T> p, const Segment<T> &s) {
    return sqrt(norm(p, s));
}

// 線分 a, b のユークリッド距離。計算量 O(1)。
template <cplib::FloatingPoint T> T distance(const Segment<T> &a, const Segment<T> &b) {
    return sqrt(norm(a, b));
}

// p のマンハッタンノルム。計算量 O(1)。
// 整数幅: B=2M。M=10^9なら32ビット、int可。
template <class T> T manhattan(Point<T> p) {
    return abs(p.x) + abs(p.y);
}

// 2点 a, b のマンハッタン距離。計算量 O(1)。
// 整数幅: B=4M。M=10^9なら33ビット、long long（int不可）。
template <class T> T manhattan(Point<T> a, Point<T> b) {
    return manhattan(a - b);
}
}
