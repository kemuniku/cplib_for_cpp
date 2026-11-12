#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 128ビット積で整数偏角を比較する。正のx軸から反時計回り。計算量 O(1)。
inline Int argcmp(std::pair<Int, Int> l, std::pair<Int, Int> r) {
    bool a = l.second > 0 || (l.second == 0 && l.first > 0),
         b = r.second > 0 || (r.second == 0 && r.first > 0);
    if (a != b)
        return b < a ? -1 : 1;
    __int128_t x = __int128_t(l.second) * r.first, y = __int128_t(l.first) * r.second;
    return x < y ? -1 : x > y ? 1 : 0;
}

// N個の非零ベクトルを正のx軸から反時計回りに安定ソートする。計算量 O(N log N)。
inline void argsort(std::vector<std::pair<Int, Int>> &x) {
    std::stable_sort(x.begin(), x.end(), [](auto l, auto r) { return argcmp(l, r) < 0; });
}

// N個の非零ベクトルを正のx軸から反時計回りに偏角ソートした列を返す。計算量 O(N log N)。
inline auto argsorted(std::vector<std::pair<Int, Int>> x) {
    argsort(x);
    return x;
}
}
