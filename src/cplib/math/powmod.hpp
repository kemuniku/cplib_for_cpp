#pragma once
#include <cplib/math/inner_math.hpp>

namespace cplib {
// 正の法について累乗の剰余を返す。O(log(n+1))。
// 正の法mについてaのn乗の剰余を0以上m未満で返す。O(log(n + 1))。
inline Int powmod(Int a, Int n, Int m) {
    assert(m > 0);
    if (m == 1)
        return 0;
    Int result = 1;
    while (n > 0) {
        if (n & 1)
            result = mul(result, a, m);
        if (n > 1)
            a = mul(a, a, m);
        n >>= 1;
    }
    return result < 0 ? result + m : result;
}
}
