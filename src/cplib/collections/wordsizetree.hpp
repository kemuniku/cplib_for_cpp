#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 64 分岐・4段の固定容量集合。更新と前後検索は O(1)。
class WordsizeTree {
    UInt a0_ = 0;
    std::array<UInt, 64> a1_{};
    std::array<UInt, 4096> a2_{};
    std::array<UInt, 262144> a3_{};
    static constexpr Int capacity = Int(1) << 24;

    static UInt highmask(Int y) {
        return y == 64 ? 0 : ~UInt(0) << y;
    }

    static UInt lowmask(Int y) {
        return y == 64 ? 0 : ~UInt(0) >> y;
    }

    template <class Range> void build(const Range &v) {
        assert(v.size() <= std::size_t(capacity));
        for (std::size_t block = 0; block < (v.size() + 63) / 64; ++block) {
            UInt bits = 0;
            for (std::size_t bit = 0; bit < std::min(std::size_t(64), v.size() - block * 64); ++bit)
                bits |= UInt(bool(v[block * 64 + bit])) << bit;
            a3_[block] = bits;
        }
        for (std::size_t i = 0; i < (v.size() + 63) / 64; ++i)
            if (a3_[i])
                a2_[i >> 6] |= UInt(1) << (i & 63);
        for (std::size_t i = 0; i < (v.size() + 4095) / 4096; ++i)
            if (a2_[i])
                a1_[i >> 6] |= UInt(1) << (i & 63);
        for (std::size_t i = 0; i < (v.size() + 262143) / 262144; ++i)
            if (a1_[i])
                a0_ |= UInt(1) << (i & 63);
    }

public:
    WordsizeTree() = default;

    explicit WordsizeTree(const std::vector<bool> &v) {
        build(v);
    }

    explicit WordsizeTree(std::span<const bool> v) {
        build(v);
    }

    void incl(Int x) {
        assert(0 <= x && x < capacity);
        a3_[x >> 6] |= UInt(1) << (x & 63);
        x >>= 6;
        a2_[x >> 6] |= UInt(1) << (x & 63);
        x >>= 6;
        a1_[x >> 6] |= UInt(1) << (x & 63);
        x >>= 6;
        a0_ |= UInt(1) << x;
    }

    bool operator[](Int x) const {
        assert(0 <= x && x < capacity);
        return (a3_[x >> 6] >> (x & 63)) & 1;
    }

    void excl(Int x) {
        assert(0 <= x && x < capacity);
        if (!(*this)[x])
            return;
        a3_[x >> 6] &= ~(UInt(1) << (x & 63));
        x >>= 6;
        if (a3_[x])
            return;
        a2_[x >> 6] &= ~(UInt(1) << (x & 63));
        x >>= 6;
        if (a2_[x])
            return;
        a1_[x >> 6] &= ~(UInt(1) << (x & 63));
        x >>= 6;
        if (a1_[x])
            return;
        a0_ &= ~(UInt(1) << x);
    }

    // x以上の最小の要素を返し、存在しなければ-1を返します。 時間計算量: O(1)、追加空間: O(1)。
    Int ge(Int x) const {
        assert(0 <= x && x < capacity);
        Int y = x & 63;
        x >>= 6;
        UInt t = a3_[x] & highmask(y);
        if (t)
            return (x << 6) | std::countr_zero(t);
        y = (x & 63) + 1;
        x >>= 6;
        t = a2_[x] & highmask(y);
        if (t) {
            x = (x << 6) | std::countr_zero(t);
            return (x << 6) | std::countr_zero(a3_[x]);
        }
        y = (x & 63) + 1;
        x >>= 6;
        t = a1_[x] & highmask(y);
        if (t) {
            x = (x << 6) | std::countr_zero(t);
            x = (x << 6) | std::countr_zero(a2_[x]);
            return (x << 6) | std::countr_zero(a3_[x]);
        }
        y = (x & 63) + 1;
        t = a0_ & highmask(y);
        if (!t)
            return -1;
        x = std::countr_zero(t);
        x = (x << 6) | std::countr_zero(a1_[x]);
        x = (x << 6) | std::countr_zero(a2_[x]);
        return (x << 6) | std::countr_zero(a3_[x]);
    }

    // x以下の最大の要素を返し、存在しなければ-1を返します。 時間計算量: O(1)、追加空間: O(1)。
    Int le(Int x) const {
        assert(0 <= x && x < capacity);
        auto last = [](UInt t) { return 63 - std::countl_zero(t); };
        Int y = 63 - (x & 63);
        x >>= 6;
        UInt t = a3_[x] & lowmask(y);
        if (t)
            return (x << 6) | last(t);
        y = 64 - (x & 63);
        x >>= 6;
        t = a2_[x] & lowmask(y);
        if (t) {
            x = (x << 6) | last(t);
            return (x << 6) | last(a3_[x]);
        }
        y = 64 - (x & 63);
        x >>= 6;
        t = a1_[x] & lowmask(y);
        if (t) {
            x = (x << 6) | last(t);
            x = (x << 6) | last(a2_[x]);
            return (x << 6) | last(a3_[x]);
        }
        y = 64 - (x & 63);
        t = a0_ & lowmask(y);
        if (!t)
            return -1;
        x = last(t);
        x = (x << 6) | last(a1_[x]);
        x = (x << 6) | last(a2_[x]);
        return (x << 6) | last(a3_[x]);
    }
};

inline WordsizeTree initWordsizeTree() {
    return {};
}

inline WordsizeTree initWordsizeTree(const std::vector<bool> &v) {
    return WordsizeTree(v);
}

inline WordsizeTree initWordsizeTree(std::span<const bool> v) {
    return WordsizeTree(v);
}
}
