#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 先頭から比較して共通接頭辞長を求める。O(min(|s|,|t|))。
template <class S, class T> Int lcp_naive(const S &s, const T &t) {
    std::size_t i = 0;
    while (i < std::min(s.size(), t.size()) && s[i] == t[i])
        ++i;
    return i;
}
}
