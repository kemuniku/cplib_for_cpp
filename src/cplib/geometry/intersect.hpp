#pragma once
#include <cplib/geometry/ccw.hpp>
#include <cplib/geometry/angle.hpp>

namespace cplib {
// 計算量 O(1)。線分a,bの交差を判定する。strict=trueで端点のみの交差を除く設定にする。
// strict=trueの共線な包含では、元実装由来で引数順に結果が変わる場合がある。
template <class T> bool intersect(const Segment<T> &a, const Segment<T> &b, bool strict = false) {
    if (strict) {
        if (ccw(a, b.s, true) == ON_SEGMENT)
            return online(a, b.t);
        if (ccw(a, b.t, true) == ON_SEGMENT)
            return online(a, b.s);
        if (ccw(b, b.s, true) == ON_SEGMENT)
            return online(b, b.t);
        if (ccw(b, b.t, true) == ON_SEGMENT)
            return online(b, b.s);
        return ccw(a, b.s) * ccw(a, b.t) < 0 && ccw(b, a.s) * ccw(b, a.t) < 0;
    }
    return ccw(a, b.s) * ccw(a, b.t) <= 0 && ccw(b, a.s) * ccw(b, a.t) <= 0;
}

// 直線 a, b が交わるかどうかを判定する。計算量 O(1)。
template <class T> bool intersect(const Line<T> &a, const Line<T> &b) {
    return !is_parallel(a, b) || online(a, b.s);
}

// 2直線 a, b の交点を返す。一致する場合は a.s を返す。計算量 O(1)。
template <class T>
    requires(!std::is_integral_v<T>)
Point<T> cross_point(const Line<T> &a, const Line<T> &b) {
    assert(intersect(a, b));
    if (is_parallel(a, b))
        return a.s;
    T d1 = cross(toVector(a), toVector(b)), d2 = cross(toVector(a), a.t - b.s);
    return b.s + toVector(b) * (d2 / d1);
}

// 2線分を延長した直線の交点を返す。計算量 O(1)。
template <class T>
    requires(!std::is_integral_v<T>)
auto cross_point(const Segment<T> &a, const Segment<T> &b) {
    return cross_point(toLine(a), toLine(b));
}

// 直線と線分を延長した直線が交わるか判定する。計算量 O(1)。
template <class T> bool intersect(const Line<T> &a, const Segment<T> &b) {
    return intersect(a, toLine(b));
}

// 線分を延長した直線と直線が交わるか判定する。計算量 O(1)。
template <class T> bool intersect(const Segment<T> &a, const Line<T> &b) {
    return intersect(toLine(a), b);
}

// 直線と線分を延長した直線の交点を返す。計算量 O(1)。
template <class T>
    requires(!std::is_integral_v<T>)
auto cross_point(const Line<T> &a, const Segment<T> &b) {
    return cross_point(a, toLine(b));
}

// 線分を延長した直線と直線の交点を返す。計算量 O(1)。
template <class T>
    requires(!std::is_integral_v<T>)
auto cross_point(const Segment<T> &a, const Line<T> &b) {
    return cross_point(toLine(a), b);
}

}
