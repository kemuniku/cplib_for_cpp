#pragma once
#include <cplib/collections/private/scalar_bitset.hpp>

namespace cplib {
using BitSet = ScalarBitSet<-1>;

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
inline BitSet initBitSetFromIndexes(const std::vector<Int> &indexes, Int n) {
    BitSet output_value(n);
    for (Int i : indexes)
        output_value.set(i, true);
    return output_value;
}
}
