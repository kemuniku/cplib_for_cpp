#pragma once
#include <cplib/common.hpp>
#include <concepts>

namespace cplib {
namespace detail {
// 桁あふれを考慮した64ビットMontgomery積。O(1)。
inline UInt prime_mont_mul(UInt a, UInt b, UInt n, UInt inverse) {
    __uint128_t t = static_cast<__uint128_t>(a) * b;
    UInt q = static_cast<UInt>(t) * inverse;
    __uint128_t sum = t + static_cast<__uint128_t>(q) * n;
    UInt r = static_cast<UInt>(sum >> 64);
    return sum < t || r >= n ? r - n : r;
}
}

// 64ビット以下の整数を決定的Miller–Rabinで判定する。O(log n)。
template <std::integral T> bool isprime(T value) {
    if (value == 2)
        return true;
    if (value < 2 || !(value & 1))
        return false;
    UInt n = value, d = n - 1;
    int s = std::countr_zero(d);
    d >>= s;
    UInt inverse = n;
    for (int i = 0; i < 6; ++i)
        inverse *= 2 - n * inverse;
    inverse = 0 - inverse;
    UInt r2 = static_cast<UInt>((-static_cast<__uint128_t>(n)) % n);
    auto mul = [&](UInt a, UInt b) { return detail::prime_mont_mul(a, b, n, inverse); };
    UInt one = mul(1, r2), minus_one = n - one;
    for (UInt a : std::array<UInt, 7>{2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {
        if (a % n == 0)
            continue;
        UInt base = mul(a % n, r2), t = one, e = d;
        while (e) {
            if (e & 1)
                t = mul(t, base);
            if (e > 1)
                base = mul(base, base);
            e >>= 1;
        }
        if (t == one || t == minus_one)
            continue;
        bool ok = false;
        for (int i = 0; i < s - 1; ++i) {
            t = mul(t, t);
            if (t == minus_one) {
                ok = true;
                break;
            }
        }
        if (!ok)
            return false;
    }
    return true;
}
}
