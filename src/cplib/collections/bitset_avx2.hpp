#pragma once
#include <cplib/collections/private/simd_bitset.hpp>

namespace cplib {
using BitSetAvx2 = SimdBitSet<2>;
}

namespace cplib {
using BitSet = BitSetAvx2;

inline BitSet initBitSet(Int n) {
    return BitSet(n);
}

inline BitSet initBitSet(const std::vector<bool> &v, Int n) {
    return BitSet(v, n);
}

inline BitSet initBitSet(const std::vector<bool> &v) {
    return BitSet(v);
}

inline BitSet initBitSetFromIndexes(const std::vector<Int> &indices, Int n) {
    BitSet out(n);
    for (Int i : indices)
        out.set(i, true);
    return out;
}

inline BitSet initBitSetFromString(std::string_view s, char match, Int n) {
    BitSet out(n);
    out.fromString(s, match);
    return out;
}

inline BitSet initBitSetFromString(std::string_view s, char match) {
    return initBitSetFromString(s, match, s.size());
}

inline BitSet initBitSetFromString(std::string_view s, std::string_view match, Int n) {
    BitSet out(n);
    out.fromString(s, match);
    return out;
}

inline BitSet initBitSetFromString(std::string_view s, std::string_view match) {
    return initBitSetFromString(s, match, s.size());
}
}
