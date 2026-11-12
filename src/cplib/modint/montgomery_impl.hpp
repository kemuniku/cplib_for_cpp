#pragma once
#include <cplib/modint/private/modint_ops.hpp>
#include <functional>

namespace cplib {
// Montgomery還元用の逆数を求める。計算量 O(1)。
constexpr std::uint32_t get_r(std::uint32_t m) {
    std::uint32_t r = m;
    for (int i = 0; i < 4; ++i)
        r *= 2u - m * r;
    return r;
}

// Montgomery変換用の 2^64 の剰余を求める。計算量 O(1)。
constexpr std::uint32_t get_n2(std::uint32_t m) {
    return std::uint32_t((~UInt(m - 1u)) % m);
}

struct MontgomeryParam {
    std::uint32_t M, r, n2;
};

template <std::uint32_t M, bool Dynamic = false>
class BasicMontgomeryModint : public detail::ModintOps<BasicMontgomeryModint<M, Dynamic>> {
    std::uint32_t a = 0;
    inline static MontgomeryParam dynamicParam{};

    static constexpr std::uint32_t reduce(UInt b) {
        auto p = get_param();
        return std::uint32_t((b + UInt(std::uint32_t(b) * (0u - p.r)) * p.M) >> 32);
    }

public:
    using montgomery_tag = void;
    using detail::ModintOps<BasicMontgomeryModint<M, Dynamic>>::operator/=;

    // 定数で割るときは逆元をコンパイル時に計算して乗算する。O(1)。
    template <Int Divisor>
    constexpr BasicMontgomeryModint &divide_constant()
        requires(!Dynamic)
    {
        constexpr auto inverse = BasicMontgomeryModint(Divisor).inv();
        return *this *= inverse;
    }

    // 還元用パラメータを返す。計算量 O(1)。
    static constexpr MontgomeryParam get_param() {
        if constexpr (Dynamic)
            return dynamicParam;
        else {
            static_assert(M < (1u << 30) && (M & 1u));
            return {M, get_r(M), get_n2(M)};
        }
    }

    // 符号なしの法を返す。計算量 O(1)。
    static constexpr std::uint32_t umod() {
        return get_param().M;
    }

    // 符号付きの法を返す。計算量 O(1)。
    static constexpr std::int32_t mod() {
        return std::int32_t(umod());
    }

    // 動的法を設定する。計算量 O(1)。
    template <std::integral I>
    static void setMod(I m)
        requires Dynamic
    {
        auto value = std::uint32_t(m);
        assert(m > 0 && UInt(m) < (UInt(1) << 30) && (value & 1u));
        dynamicParam = {value, get_r(value), get_n2(value)};
        assert(dynamicParam.r * value == 1u);
    }

    // 零を構築する。計算量 O(1)。
    constexpr BasicMontgomeryModint() = default;

    // 整数の剰余を構築する。計算量 O(1)。
    template <std::integral I> constexpr BasicMontgomeryModint(I value) {
        auto p = get_param();
        a = reduce(detail::normalizeModint(value, p.M) * p.n2);
    }

    // [0,2M)の冗長表現で加減乗算する。値取得時のみ[0,M)へ正規化。
    // 剰余を加算する。計算量 O(1)。
    constexpr BasicMontgomeryModint &operator+=(BasicMontgomeryModint b) {
        a += b.a - umod() * 2u;
        if (std::bit_cast<std::int32_t>(a) < 0)
            a += umod() * 2u;
        return *this;
    }

    // 剰余を減算する。計算量 O(1)。
    constexpr BasicMontgomeryModint &operator-=(BasicMontgomeryModint b) {
        a -= b.a;
        if (std::bit_cast<std::int32_t>(a) < 0)
            a += umod() * 2u;
        return *this;
    }

    // 剰余を乗算する。計算量 O(1)。
    constexpr BasicMontgomeryModint &operator*=(BasicMontgomeryModint b) {
        a = reduce(UInt(a) * b.a);
        return *this;
    }

    // 加法逆元を返す。計算量 O(1)。
    constexpr BasicMontgomeryModint operator-() const {
        BasicMontgomeryModint result = 0;
        result -= *this;
        return result;
    }

    // 正規化された剰余を返す。計算量 O(1)。
    constexpr Int val() const {
        Int result = reduce(a);
        if (UInt(result) >= umod())
            result -= umod();
        return result;
    }

    // 冗長な内部表現を正規化して剰余の等値を判定する。O(1)。
    // 剰余が等しいか判定する。計算量 O(1)。
    constexpr bool operator==(const BasicMontgomeryModint &b) const {
        auto m = umod();
        return (a >= m ? a - m : a) == (b.a >= m ? b.a - m : b.a);
    }
};

template <std::uint32_t M> using StaticMontgomeryModint = BasicMontgomeryModint<M, false>;
template <std::uint32_t Id> using DynamicMontgomeryModint = BasicMontgomeryModint<Id, true>;
template <class T>
concept MontgomeryModint = requires { typename T::montgomery_tag; };
}

namespace std {
// 同じ剰余が同じハッシュ値になるように計算する。O(1)。
template <uint32_t M, bool D> struct hash<cplib::BasicMontgomeryModint<M, D>> {
    // 剰余のハッシュ値を返す。計算量 O(1)。
    size_t operator()(cplib::BasicMontgomeryModint<M, D> x) const {
        return hash<cplib::Int>{}(x.val());
    }
};
}

#define declarStaticMontgomeryModint(NAME, M) using NAME = ::cplib::StaticMontgomeryModint<M>
#define declarDynamicMontgomeryModint(NAME, ID) using NAME = ::cplib::DynamicMontgomeryModint<ID>
