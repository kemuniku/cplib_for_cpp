#pragma once
#include <cplib/common.hpp>

namespace cplib {
using float128 = __float128;

// 128ビット浮動小数点数に変換する。O(1)。
template <class T> constexpr float128 to_float128(T x) {
    return static_cast<float128>(x);
}

// double に変換する。計算量 O(1)。
inline double to_float(float128 x) {
    return static_cast<double>(x);
}

// 絶対値を返す。計算量 O(1)。
inline float128 abs(float128 x) {
    return x >= 0 ? x : -x;
}

// 大小比較し、小・等・大に対して -1・0・1 を返す。計算量 O(1)。
inline int cmp(float128 x, float128 y) {
    return x < y ? -1 : x == y ? 0 : 1;
}
}
