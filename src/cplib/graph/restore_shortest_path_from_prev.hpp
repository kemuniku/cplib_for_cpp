#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 前駆配列から経路を復元する。O(経路長)。到達不能なgoalも単独で返す。
inline std::vector<Int> restore_shortest_path_from_prev(std::span<const Int> prev, Int goal) {
    std::vector<Int> result;
    for (Int i = goal; i != -1; i = prev[i])
        result.push_back(i);
    std::reverse(result.begin(), result.end());
    return result;
}
}
