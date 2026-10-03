#pragma once
#include <cplib/common.hpp>

namespace cplib {
// ax+by=gcd(a,b) の係数を元実装と同じ互除法で返す。O(log(max(|a|,|b|)))。
inline std::pair<Int, Int> ext_gcd(Int aa, Int bb) {
    bool ia = aa < 0, ib = bb < 0;
    __int128 a = aa, b = bb;
    if (ia)
        a = -a;
    if (ib)
        b = -b;
    bool sw = a < b;
    if (sw)
        std::swap(a, b);
    std::vector<std::pair<__int128, __int128>> line{{a, b}};
    while (a > 0 && b > 0) {
        if (line.size() & 1)
            a %= b;
        else
            b %= a;
        line.emplace_back(a, b);
    }
    __int128 x = line.size() & 1, y = 1 ^ (line.size() & 1);
    for (std::size_t i = line.size() - 1; i > 0;) {
        auto [u, v] = line[--i];
        if (u < v)
            x -= v / u * y;
        else
            y -= u / v * x;
    }
    if (sw)
        std::swap(x, y);
    if (ia)
        x = -x;
    if (ib)
        y = -y;
    return {static_cast<Int>(x), static_cast<Int>(y)};
}
}
