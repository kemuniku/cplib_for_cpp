#pragma once
#include <cplib/collections/private/scalar_bitset.hpp>

namespace cplib {
template <Int N> using BitSet = ScalarBitSet<N>;

// 長さNの全ビット0の列を構築する。O(1+N/64)。
template <Int N> BitSet<N> initBitSet() {
    return {};
}

// 真偽値配列から長さNのビット列を構築する。O(1+N/64+v.size())。
template <Int N> BitSet<N> initBitSet(const std::vector<bool> &v) {
    return BitSet<N>(v);
}

// 指定した添字のみ1の長さNの列を構築する。添字数kとしてO(1+N/64+k)。
template <Int N> BitSet<N> initBitSetFromIndexes(const std::vector<Int> &indices) {
    BitSet<N> out;
    for (Int i : indices)
        out.set(i, true);
    return out;
}
}
