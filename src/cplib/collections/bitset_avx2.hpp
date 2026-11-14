#pragma once
#include <cplib/collections/private/simd_bitset.hpp>

namespace cplib {
using BitSetAvx2 = SimdBitSet<2>;
}

namespace cplib::avx2 {
using BitSet = BitSetAvx2;

// 長さnの全ビット0の列を構築する。O(1+n/64)。
inline BitSet initBitSet(Int n) {
    return BitSet(n);
}

// 真偽値配列からビット列を構築する。結果の長さをnとしてO(1+n/64+v.size())。
inline BitSet initBitSet(const std::vector<bool> &v, Int n) {
    return BitSet(v, n);
}

// 真偽値配列からビット列を構築する。結果の長さをnとしてO(1+n/64+v.size())。
inline BitSet initBitSet(const std::vector<bool> &v) {
    return BitSet(v);
}

// 指定した添字のみ1の長さnの列を構築する。添字数をkとしてO(1+n/64+k)。
inline BitSet initBitSetFromIndexes(const std::vector<Int> &indices, Int n) {
    BitSet output_value(n);
    for (Int i : indices)
        output_value.set(i, true);
    return output_value;
}

// 文字列sでmatchと一致する文字の位置を1にする。結果の長さをnとしてO(1+n/64+s.size())。
inline BitSet initBitSetFromString(std::string_view s, char match, Int n) {
    BitSet output_value(n);
    output_value.fromString(s, match);
    return output_value;
}

// 文字列sでmatchと一致する文字の位置を1にする。結果の長さをnとしてO(1+n/64+s.size())。
inline BitSet initBitSetFromString(std::string_view s, char match) {
    return initBitSetFromString(s, match, s.size());
}

// 文字列sでmatchの同じ位置の文字と一致する文字の位置を1にする。結果の長さをnとしてO(1+n/64+s.size())。
inline BitSet initBitSetFromString(std::string_view s, std::string_view match, Int n) {
    BitSet output_value(n);
    output_value.fromString(s, match);
    return output_value;
}

// 文字列sでmatchの同じ位置の文字と一致する文字の位置を1にする。結果の長さをnとしてO(1+n/64+s.size())。
inline BitSet initBitSetFromString(std::string_view s, std::string_view match) {
    return initBitSetFromString(s, match, s.size());
}
}
