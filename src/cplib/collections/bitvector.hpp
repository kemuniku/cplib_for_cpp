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
    Int length_;
    std::vector<UInt> bits;
    std::vector<Int> csum;

public:
    // lengthビットの零ビット列を構築する。O(length/64+1)。
    explicit BitVector(Int length = 0)
        : length_(length), bits((length + 63) / 64 + 1), csum((length + 63) / 64 + 1) {
    }

    // パディングを含まない論理ビット数。O(1)。
    Int size() const {
        return length_;
    }

    // idxビットを1にする。O(1)。
    // buildする前にだけ呼ぶ
    void set(Int idx) {
        bits[idx >> 6] |= UInt(1) << (idx & 63);
    }

    // 計算量はO(1)。
    // 長さ外のビットは呼出側で0にする。設定後にbuildする。
    void setWord(Int idx, UInt value) {
        bits[idx] = value;
    }

    // rank用の累積和を構築する。O(length/64+1)。
    void build() {
        for (std::size_t i = 0; i + 1 < bits.size(); ++i)
            csum[i + 1] = csum[i] + detail::bitvector_popcount(bits[i]);
    }

    // idx番目のビットを返す。O(1)。
    bool access(Int idx) const {
        return (bits[idx >> 6] >> (idx & 63)) & 1;
    }

    // idx番目のビットを返す。O(1)。
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

// lengthビットの零ビット列を生成する。O(length/64+1)。
inline auto newBitVector(Int length) {
    return BitVector(length);
}

// ビルド前にi番目のビットを1にする。O(1)。
inline void set(BitVector &b, Int i) {
    b.set(i);
}

// ビルド前にi番目のワードを設定する。O(1)。
inline void setWord(BitVector &b, Int i, UInt value) {
    b.setWord(i, value);
}

// rank用の累積和を構築する。ビット数をnとしてO(n/64+1)。
inline void build(BitVector &b) {
    b.build();
}

// i番目のビットを返す。O(1)。
inline bool access(const BitVector &b, Int i) {
    return b[i];
}

// [0,i)に含まれる1の個数を返す。O(1)。
inline Int rank(const BitVector &b, Int i) {
    return b.rank(i);
}
}
