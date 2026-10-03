#pragma once
#include <cplib/geometry/intersect.hpp>

namespace cplib {
// 点 a, b のユークリッド距離の2乗
template <class T> T norm(Point<T> a, Point<T> b) {
    return norm(a - b);
}

// 点 p と直線 l のユークリッド距離の2乗
template <class T> T norm(Point<T> p, const Line<T> &l) {
    T area = cross(vector(l), p - l.s);
    return area * area / norm(vector(l));
}

// 点 p と線分 s のユークリッド距離の2乗
template <class T> T norm(Point<T> p, const Segment<T> &s) {
    if (geometry_lt(dot(vector(s), p - s.s), 0))
        return norm(p - s.s);
    if (geometry_lt(dot(-vector(s), p - s.t), 0))
        return norm(p - s.t);
    return norm(p, toLine(s));
}

// 線分 a, b のユークリッド距離の2乗
template <class T> T norm(const Segment<T> &a, const Segment<T> &b) {
    if (intersect(a, b))
        return T(0);
    return std::min({norm(a.s, b), norm(a.t, b), norm(b.s, a), norm(b.t, a)});
}

// 点 a, b のユークリッド距離
template <std::floating_point T> double distance(Point<T> a, Point<T> b) {
    return std::sqrt(norm(a, b));
}

// 点 p と直線 l のユークリッド距離
template <std::floating_point T> double distance(Point<T> p, const Line<T> &l) {
    return std::sqrt(norm(p, l));
}

// 点 p と線分 s のユークリッド距離
template <std::floating_point T> double distance(Point<T> p, const Segment<T> &s) {
    return std::sqrt(norm(p, s));
}

// 線分 a, b のユークリッド距離
template <std::floating_point T> double distance(const Segment<T> &a, const Segment<T> &b) {
    return std::sqrt(norm(a, b));
}

// p のマンハッタンノルム
template <class T> T manhattan(Point<T> p) {
    using std::abs;
    return abs(p.x) + abs(p.y);
}

// 2点 a, b のマンハッタン距離
template <class T> T manhattan(Point<T> a, Point<T> b) {
    return manhattan(a - b);
}
}
