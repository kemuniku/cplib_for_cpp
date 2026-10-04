#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 非負整数の平方根を切り捨てて返す。Newton法、O(log n)。
inline Int isqrt(Int n) {
    assert(n >= 0);
    Int x = n, y = x / 2 + x % 2;
    while (y < x) {
        x = y;
        y = static_cast<Int>((static_cast<__int128>(x) + n / x) / 2);
    }
    return x;
}
}
