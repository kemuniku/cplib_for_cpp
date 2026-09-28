#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 非負整数の最小欠損値を返す。O(n)時間・領域。
template<class Range> Int mex_naive(const Range& xs) {
    std::vector<bool> used(xs.size());
    for (auto x:xs) if (x>=0 && static_cast<std::size_t>(x)<used.size()) used[x]=true;
    for (std::size_t i=0;i<used.size();++i) if (!used[i]) return i;
    return used.size();
}
}
