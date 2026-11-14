#pragma once
#include <cplib/geometry/base.hpp>

namespace cplib {
// 直線 l への点 p の射影。計算量 O(1)。整数座標では途中の除算が整数除算になる。
template <class T> Point<T> projection(const Line<T> &l, Point<T> p) {
    T t = dot(p - l.s, toVector(l)) / norm(toVector(l));
    return l.s + toVector(l) * t;
}

// 点 p の直線 l に対する反射。計算量 O(1)。整数座標では射影の途中が整数除算になる。
template <class T> Point<T> reflection(const Line<T> &l, Point<T> p) {
    return p + (projection(l, p) - p) * 2;
}

// 線分を延長した直線への点pの射影を返す。計算量 O(1)。整数座標では途中が整数除算になる。
template <class T> auto projection(const Segment<T> &l, Point<T> p) {
    return projection(toLine(l), p);
}

// 線分を延長した直線に対する点pの反射を返す。計算量 O(1)。整数座標では途中が整数除算になる。
template <class T> auto reflection(const Segment<T> &l, Point<T> p) {
    return reflection(toLine(l), p);
}
}
