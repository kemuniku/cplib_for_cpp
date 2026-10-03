#pragma once
#include <cplib/geometry/base.hpp>

namespace cplib {
// 直線 l への点 p の射影
template <class T> Point<T> projection(const Line<T> &l, Point<T> p) {
    T t = dot(p - l.s, vector(l)) / norm(vector(l));
    return l.s + vector(l) * t;
}

// 点 p の直線 l に対する反射
template <class T> Point<T> reflection(const Line<T> &l, Point<T> p) {
    return p + (projection(l, p) - p) * 2;
}

template <class T> auto projection(const Segment<T> &l, Point<T> p) {
    return projection(toLine(l), p);
}

template <class T> auto reflection(const Segment<T> &l, Point<T> p) {
    return reflection(toLine(l), p);
}
}
