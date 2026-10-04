#pragma once
#include <cplib/common.hpp>

namespace cplib {
class XorBasis {
    // 基底を簡約し、新しい独立要素を追加する。O(基底数)。
    void insert_reduced(Int e) {
        for (Int b : basis)
            if ((e ^ b) < e)
                e ^= b;
        if (e) {
            for (Int &b : basis)
                if ((e ^ b) < b)
                    b ^= e;
            basis.push_back(e);
        }
    }

public:
    std::vector<Int> basis;
    XorBasis() = default;

    // 簡約基底を構築する。O(n log V)。
    explicit XorBasis(std::span<const Int> a) {
        for (Int e : a)
            insert_reduced(e);
        std::sort(basis.begin(), basis.end());
    }

    // 値を追加する。O(log V log log V)。
    void incl(Int x) {
        insert_reduced(x);
        std::sort(basis.begin(), basis.end());
    }

    Int len_basis() const {
        return basis.size();
    }

    // 基底のXOR和でxが作成可能か？
    bool can_make(Int x) const {
        for (Int b : basis)
            if ((x ^ b) < x)
                x ^= b;
        return x == 0;
    }

    // 作成可能な値のk番目を返す。O(log V)。
    // 作成可能な値の中でk番目に小さいものを求めます。
    // 存在しない場合は、-1を返します。
    Int kth_smallest(Int k) const {
        if (k < 0 || (basis.size() < 64 && UInt(k) >= (UInt(1) << basis.size())))
            return -1;
        Int result = 0;
        for (std::size_t i = 0; i < basis.size(); ++i)
            if (UInt(k) & (UInt(1) << i))
                result ^= basis[i];
        return result;
    }

    // 作成可能な値の中でx未満の内最大のもの
    // 存在しない場合は-1
    Int lt(Int x) const {
        if (x == 0)
            return -1;
        Int r = 0;
        for (auto it = basis.rbegin(); it != basis.rend(); ++it)
            if ((r ^ *it) < x)
                r ^= *it;
        return r;
    }

    // 作成可能な値の中でx以下の内最大のもの
    Int le(Int x) const {
        Int r = 0;
        for (auto it = basis.rbegin(); it != basis.rend(); ++it)
            if ((r ^ *it) <= x)
                r ^= *it;
        return r;
    }

    // 作成可能な値の中でxが小さい方から何番目かを返す。結果の各ビットは使用する基底を表す。
    Int index(Int x) const {
        UInt r = 0;
        for (std::size_t i = 0; i < basis.size(); ++i)
            if ((x ^ basis[i]) < x) {
                x ^= basis[i];
                r += UInt(1) << i;
            }
        return static_cast<Int>(r);
    }

    // xとのXORが最小になる作成可能値を返す。O(log V)。
    Int xor_min(Int x) const {
        Int r = x;
        for (auto it = basis.rbegin(); it != basis.rend(); ++it)
            if ((r ^ *it) < r)
                r ^= *it;
        return r ^ x;
    }

    // 作成可能な値の中でxとxorを取ったときにk番目に小さくなるもの
    // 存在しないならば-1を返す。
    Int xor_kth(Int x, Int k) const {
        Int kth = kth_smallest(k);
        if (kth == -1)
            return -1;
        Int v = xor_min(x) ^ x;
        for (auto it = basis.rbegin(); it != basis.rend(); ++it)
            if ((v ^ *it) < v)
                v ^= *it;
        return v ^ kth ^ x;
    }
};

inline XorBasis initXorBasis(std::span<const Int> a) {
    return XorBasis(a);
}

inline XorBasis initXorBasis() {
    return {};
}

inline void incl(XorBasis &s, Int x) {
    s.incl(x);
}

inline Int len_basis(const XorBasis &s) {
    return s.len_basis();
}

inline bool can_make(const XorBasis &s, Int x) {
    return s.can_make(x);
}

inline Int kth_smallest(const XorBasis &s, Int k) {
    return s.kth_smallest(k);
}

inline Int lt(const XorBasis &s, Int x) {
    return s.lt(x);
}

inline Int le(const XorBasis &s, Int x) {
    return s.le(x);
}

inline Int index(const XorBasis &s, Int x) {
    return s.index(x);
}

inline Int xor_min(const XorBasis &s, Int x) {
    return s.xor_min(x);
}

inline Int xor_kth(const XorBasis &s, Int x, Int k) {
    return s.xor_kth(x, k);
}
}
