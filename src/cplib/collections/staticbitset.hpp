#pragma once
#include <cplib/collections/private/scalar_bitset.hpp>

namespace cplib {
template <Int N> using BitSet = ScalarBitSet<N>;

template <Int N> BitSet<N> initBitSet() {
    return {};
}

template <Int N> BitSet<N> initBitSet(const std::vector<bool> &v) {
    return BitSet<N>(v);
}

template <Int N> BitSet<N> initBitSetFromIndexes(const std::vector<Int> &indices) {
    BitSet<N> out;
    for (Int i : indices)
        out.set(i, true);
    return out;
}
}
