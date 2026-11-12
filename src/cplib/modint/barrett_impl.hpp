#pragma once
#include <cplib/modint/private/modint_ops.hpp>
#include <functional>

namespace cplib {
// Barrett還元用の逆数を返す。計算量 O(1)。
constexpr UInt get_im(std::uint32_t m) {
    return std::numeric_limits<UInt>::max() / m + 1;
}

// 64ビット値の積の上位64ビットを返す。計算量 O(1)。
constexpr UInt calc_mul(UInt a, UInt b) {
    return UInt((__uint128_t(a) * b) >> 64);
}

struct BarrettParam {
    std::uint32_t M;
    UInt im;
};

template <std::uint32_t M, bool Dynamic = false>
class BasicBarrettModint : public detail::ModintOps<BasicBarrettModint<M, Dynamic>> {
    std::uint32_t a = 0;
    inline static BarrettParam dynamicParam{};

public:
    using barrett_tag = void;
    using detail::ModintOps<BasicBarrettModint<M, Dynamic>>::operator/=;

    // 定数で割るときは逆元をコンパイル時に計算して乗算する。O(1)。
    template <Int Divisor>
    constexpr BasicBarrettModint &divide_constant()
        requires(!Dynamic)
    {
        constexpr auto inverse = BasicBarrettModint(Divisor).inv();
        return *this *= inverse;
    }

    // 還元用パラメータを返す。計算量 O(1)。
    static constexpr BarrettParam get_param() {
        if constexpr (Dynamic)
            return dynamicParam;
        else {
            static_assert(M > 0 && M <= std::uint32_t(std::numeric_limits<std::int32_t>::max()));
            return {M, get_im(M)};
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
        assert(m > 0 && UInt(m) <= UInt(std::numeric_limits<std::int32_t>::max()));
        dynamicParam = {std::uint32_t(m), get_im(std::uint32_t(m))};
    }

    // 64bit積の上位を用いるBarrett還元。O(1)。
    static constexpr std::uint32_t rem(UInt value) {
        if constexpr (!Dynamic)
            return std::uint32_t(value % M);
        auto p = get_param();
        UInt x = calc_mul(value, p.im), r = value - x * p.M;
        if (p.M <= r)
            r += p.M;
        return std::uint32_t(r);
    }

    // 零を構築する。計算量 O(1)。
    constexpr BasicBarrettModint() = default;

    // 整数の剰余を構築する。計算量 O(1)。
    template <std::integral I>
    constexpr BasicBarrettModint(I value) : a(detail::normalizeModint(value, umod())) {
    }

    // 正規化された剰余を返す。計算量 O(1)。
    constexpr Int val() const {
        return a;
    }

    // 加法逆元を返す。計算量 O(1)。
    constexpr BasicBarrettModint operator-() const {
        auto result = *this;
        if (a)
            result.a = umod() - a;
        return result;
    }

    // 剰余を加算する。計算量 O(1)。
    constexpr BasicBarrettModint &operator+=(BasicBarrettModint b) {
        a += b.a;
        if (a >= umod())
            a -= umod();
        return *this;
    }

    // 剰余を減算する。計算量 O(1)。
    constexpr BasicBarrettModint &operator-=(BasicBarrettModint b) {
        a -= b.a;
        if (a >= umod())
            a += umod();
        return *this;
    }

    // 剰余を乗算する。計算量 O(1)。
    constexpr BasicBarrettModint &operator*=(BasicBarrettModint b) {
        a = rem(UInt(a) * b.a);
        return *this;
    }

    // 剰余が等しいか判定する。計算量 O(1)。
    constexpr bool operator==(const BasicBarrettModint &b) const {
        return a == b.a;
    }
};

template <std::uint32_t M> using StaticBarrettModint = BasicBarrettModint<M, false>;
template <std::uint32_t Id> using DynamicBarrettModint = BasicBarrettModint<Id, true>;
template <class T>
concept BarrettModint = requires { typename T::barrett_tag; };

// 指定型の法で剰余を求める。計算量 O(1)。
template <BarrettModint T> constexpr auto rem(UInt value) {
    return T::rem(value);
}
}

namespace std {
template <uint32_t M, bool D> struct hash<cplib::BasicBarrettModint<M, D>> {
    // 剰余のハッシュ値を返す。計算量 O(1)。
    size_t operator()(cplib::BasicBarrettModint<M, D> x) const {
        return hash<cplib::Int>{}(x.val());
    }
};
}

#define declarStaticBarrettModint(NAME, M) using NAME = ::cplib::StaticBarrettModint<M>
#define declarDynamicBarrettModint(NAME, ID) using NAME = ::cplib::DynamicBarrettModint<ID>
