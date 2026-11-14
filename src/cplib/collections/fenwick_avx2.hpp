#pragma once
#include <cplib/collections/private/fenwick_avx2_kernel.hpp>
#include <cplib/utils/backwards_index.hpp>

namespace cplib {
class FenwickTreeAvx2 {
    Int size_, height = 0;
    std::array<Int, 16> offsets{};
    std::vector<Int> data;

public:
    // 16分岐の部分累積和。O(N)構築、約16N/15個の64ビット整数。
    explicit FenwickTreeAvx2(Int n = 0) : size_(n) {
        assert(n >= 0);
        Int m = n, total = 0;
        while (true) {
            offsets[height++] = total;
            total += (m + 16) & ~Int(15);
            m >>= 4;
            if (m == 0)
                break;
        }
        data.resize(total + 15);
        Int alignment = ((UInt(0) - reinterpret_cast<std::uintptr_t>(data.data())) & 127) >> 3;
        for (Int h = 0; h < height; ++h)
            offsets[h] += alignment;
    }

    // 配列から構築する。要素の演算がO(1)ならO(n)。
    explicit FenwickTreeAvx2(std::span<const Int> values) : FenwickTreeAvx2(values.size()) {
        std::copy(values.begin(), values.end(), data.begin() + offsets[0]);
        detail::fenwick_avx2::cplib_fw16_build(data.data(), offsets.data(), height, size_);
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return size_;
    }

    // p番目の値にdeltaを加算する。O(log₁₆(n+1))。
    void add(Int p, Int delta) {
        assert(0 <= p && p < size_);
        detail::fenwick_avx2::cplib_fw16_add(data.data(), offsets.data(), height, p, delta);
    }

    // [0,r)の和をIntとしてO(log_16 n)で返します。
    Int prefix(Int r) const {
        assert(0 <= r && r <= size_);
        return r == 0 ? 0 : detail::fenwick_avx2::cplib_fw16_prefix(data.data(), offsets.data(), r);
    }

    // [l,r)の和をIntとしてO(log_16 n)で返します。
    Int get(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= size_);
        return l == r ? 0 : detail::fenwick_avx2::cplib_fw16_get(data.data(), offsets.data(), l, r);
    }

    // p番目の値を返す。O(log₁₆(n+1))。
    Int operator[](Int p) const {
        return get(p, p + 1);
    }

    // 両端を含む区間の和を返す。O(log₁₆(n+1))。
    Int operator[](ClosedSlice<Int, Int> s) const {
        return get(s.a, s.b + 1);
    }

    // p番目の値を設定する。O(log₁₆(n+1))。
    void set(Int p, Int value) {
        add(p, detail::fenwick_avx2::cplib_fw16_diff(value, get(p, p + 1)));
    }

    struct Reference {
        FenwickTreeAvx2 *owner;
        Int index;

        // 指す位置の値を取得する。O(log₁₆(n+1))。
        operator Int() const {
            return std::as_const(*owner)[index];
        }

        // 指す位置の値を設定する。O(log₁₆(n+1))。
        Reference &operator=(Int v) {
            owner->set(index, v);
            return *this;
        }

        // 別の位置の値を読み取り代入する。O(log₁₆(n+1))。
        Reference &operator=(const Reference &other) {
            return *this = Int(other);
        }
    };

    // p番目の要素を読み書きする参照を返す。O(1)。
    Reference operator[](Int p) {
        return {this, p};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
};

// 長さnの零配列から構築する。O(n)。
inline auto initFenwickTreeAvx2(Int n) {
    return FenwickTreeAvx2(n);
}

// 配列から構築する。長さをnとしてO(n)。
inline auto initFenwickTreeAvx2(std::span<const Int> v) {
    return FenwickTreeAvx2(v);
}

// p番目に加算する。O(log₁₆(n+1))。
inline void add(FenwickTreeAvx2 &t, Int p, Int delta) {
    t.add(p, delta);
}

// [0,r)の和を返す。O(log₁₆(n+1))。
inline Int prefix(const FenwickTreeAvx2 &t, Int r) {
    return t.prefix(r);
}

// [l,r)の和を返す。O(log₁₆(n+1))。
inline Int get(const FenwickTreeAvx2 &t, Int l, Int r) {
    return t.get(l, r);
}
}
