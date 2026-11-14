#pragma once
#include <cplib/collections/private/simd_bitset.hpp>

namespace cplib::static_bitset_avx2 {
template <Int N> using BitSet = SimdBitSet<2, N>;

// 長さNの全ビット0の列を構築する。O(1+N/64)。
template <Int N> BitSet<N> initBitSet() {
    return {};
}

// 真偽値配列から長さNのビット列を構築する。O(1+N/64+v.size())。
template <Int N> BitSet<N> initBitSet(const std::vector<bool> &v) {
    return BitSet<N>(v);
}

// 真偽値配列から長さNのビット列を構築する。O(1+N/64+v.size())。
template <Int N> BitSet<N> initBitSet(std::span<const bool> v) {
    return BitSet<N>(v);
}

// 指定した添字のみ1の長さNの列を構築する。添字数kとしてO(1+N/64+k)。
template <Int N> BitSet<N> initBitSetFromIndexes(const std::vector<Int> &indices) {
    BitSet<N> out;
    for (Int i : indices)
        out.set(i, true);
    return out;
}

// 文字列sでmatchと一致する文字の位置を1にする。O(1+N/64+s.size())。
template <Int N> BitSet<N> initBitSetFromString(std::string_view s, char match) {
    BitSet<N> out;
    out.fromString(s, match);
    return out;
}

// 文字列sでmatchの同じ位置の文字と一致する文字の位置を1にする。O(1+N/64+s.size())。
template <Int N> BitSet<N> initBitSetFromString(std::string_view s, std::string_view match) {
    BitSet<N> out;
    out.fromString(s, match);
    return out;
}
}
