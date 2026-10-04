#pragma once
#include <cplib/common.hpp>

namespace cplib {
namespace detail {
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("popcnt")))
#endif
inline Int bitvector_popcount(UInt x) {
    return std::popcount(x);
}
}

class BitVector {
    std::vector<UInt> bits;
    std::vector<Int> csum;

public:
    explicit BitVector(Int length = 0)
        : bits((length + 63) / 64 + 1), csum((length + 63) / 64 + 1) {
    }

    // buildする前にだけ呼ぶ
    void set(Int idx) {
        bits[idx >> 6] |= UInt(1) << (idx & 63);
    }

    // 長さ外のビットは呼出側で0にする。設定後にbuildする。
    void setWord(Int idx, UInt value) {
        bits[idx] = value;
    }

    void build() {
        for (std::size_t i = 0; i + 1 < bits.size(); ++i)
            csum[i + 1] = csum[i] + detail::bitvector_popcount(bits[i]);
    }

    bool access(Int idx) const {
        return (bits[idx >> 6] >> (idx & 63)) & 1;
    }

    bool operator[](Int idx) const {
        return access(idx);
    }

    // [0,idx)の1の個数。O(1)。
    // [0,idx)の1の個数をO(1)で返します。build後に呼んでください。
    Int rank(Int idx) const {
        return csum[idx >> 6] +
               detail::bitvector_popcount(bits[idx >> 6] & ((UInt(1) << (idx & 63)) - 1));
    }
};

inline auto newBitVector(Int length) {
    return BitVector(length);
}

inline void set(BitVector &b, Int i) {
    b.set(i);
}

inline void setWord(BitVector &b, Int i, UInt value) {
    b.setWord(i, value);
}

inline void build(BitVector &b) {
    b.build();
}

inline bool access(const BitVector &b, Int i) {
    return b[i];
}

inline Int rank(const BitVector &b, Int i) {
    return b.rank(i);
}
}
