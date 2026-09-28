#pragma once
#include <cplib/common.hpp>
namespace cplib {
// (gcd(a,b), 逆元係数)を返す。b>=1、O(log b)。
inline std::pair<Int, Int> inv_gcd(Int a, Int b) {
    assert(b >= 1);
    a %= b;
    if (a < 0) a += b;
    if (!a) return {b, 0};
    Int s = b, t = a, m0 = 0, m1 = 1;
    while (t) {
        Int u = s / t;
        s -= t * u; m0 -= m1 * u;
        std::swap(s, t); std::swap(m0, m1);
    }
    if (m0 < 0) m0 += b / s;
    return {s, m0};
}
}
