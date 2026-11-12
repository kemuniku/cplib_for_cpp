#pragma once
#include <cplib/math/inv_gcd.hpp>
#include <stdexcept>

namespace cplib {
// 一般の合同式を併合し最小非負剰余と最小公倍数を返す。解なしは(0,0)、空は(0,1)。
// 不正入力はinvalid_argument、途中の法の桁あふれはoverflow_error。時間 O(N log(max m))、領域 O(1)。
inline std::pair<Int, Int> crt(const std::vector<Int> &r, const std::vector<Int> &m) {
    if (r.size() != m.size())
        throw std::invalid_argument("剰余列と法列の長さが一致しません");
    for (Int modulus : m)
        if (modulus <= 0)
            throw std::invalid_argument("CRTの法は正の整数が必要です");
    Int residue = 0, modulus = 1;
    for (std::size_t i = 0; i < r.size(); ++i) {
        Int next = r[i] % m[i];
        if (next < 0)
            next += m[i];
        auto [g, inverse] = inv_gcd(modulus, m[i]);
        Int difference = next - residue;
        if (difference % g != 0)
            return {0, 0};
        Int factor = m[i] / g;
        if (modulus > std::numeric_limits<Int>::max() / factor)
            throw std::overflow_error("CRTの最小公倍数がIntの範囲外です");
        __int128 step = __int128(difference / g) * inverse % factor;
        if (step < 0)
            step += factor;
        residue += modulus * Int(step);
        modulus *= factor;
    }
    return {residue, modulus};
}
}
