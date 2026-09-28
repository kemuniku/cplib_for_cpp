#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 連続する同じ値を(値,個数)に圧縮する。O(n)。
template<class Range> auto run_length_encode(const Range& a) {
    using T=typename Range::value_type;
    std::vector<std::pair<T,Int>> result;
    for (const auto& x:a) {
        if (!result.empty() && result.back().first==x) ++result.back().second;
        else result.emplace_back(x,1);
    }
    return result;
}
}
