#pragma once
#include <cplib/collections/private/simd_bitset.hpp>

namespace cplib {
template <Int N> using BitSet = SimdBitSet<512, N>;

template <Int N> BitSet<N> initBitSet() {
    return {};
}

template <Int N> BitSet<N> initBitSet(const std::vector<bool> &v) {
    return BitSet<N>(v);
}

template <Int N> BitSet<N> initBitSet(std::span<const bool> v) {
    return BitSet<N>(v);
}

template <Int N> BitSet<N> initBitSetFromIndexes(const std::vector<Int> &indices) {
    BitSet<N> out;
    for (Int i : indices)
        out.set(i, true);
    return out;
}

template <Int N> BitSet<N> initBitSetFromString(std::string_view s, char match) {
    BitSet<N> out;
    out.fromString(s, match);
    return out;
}

template <Int N> BitSet<N> initBitSetFromString(std::string_view s, std::string_view match) {
    BitSet<N> out;
    out.fromString(s, match);
    return out;
}
}
