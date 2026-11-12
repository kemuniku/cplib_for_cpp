#pragma once
#include <cplib/common.hpp>
#include <ostream>

namespace cplib {
namespace detail {
constexpr Int expmod_phi(Int n) {
    Int result = n;
    for (Int p = 2; p <= n / p; ++p)
        if (n % p == 0) {
            result -= result / p;
            while (n % p == 0)
                n /= p;
        }
    if (n > 1)
        result -= result / n;
    return result;
}

template <Int P> constexpr Int expmod_getmod(Int x) {
    return x < 2 * P ? x : x % P + P;
}

template <Int P> constexpr Int expmod_mul(Int x, Int y) {
    return expmod_getmod<P>(x * y);
}

template <Int P> constexpr Int expmod_pow(Int a, Int n) {
    Int result = 1 % P;
    while (n > 0) {
        if (n & 1)
            result = expmod_mul<P>(result, a);
        if (n > 1)
            a = expmod_mul<P>(a, a);
        n >>= 1;
    }
    return result;
}
}

// phi鎖と[P,2P)の印を保持する剰余。積がIntを越えない入力を前提とする。
// 以下 L は P から 1 までオイラーのφ関数を繰り返す鎖の長さ。
template <Int P> struct expmodint {
    static_assert(P > 1);
    Int x = 0;
    expmodint<detail::expmod_phi(P)> p;
    // 零を構築する。計算量 O(L)。
    constexpr expmodint() = default;

    // 整数から剰余とφ鎖を構築する。計算量 O(L)。
    constexpr explicit expmodint(Int value) : x(detail::expmod_getmod<P>(value)), p(value) {
    }

    // 現在の値と子のφ鎖から構築する。計算量 O(L)。
    constexpr expmodint(Int value, expmodint<detail::expmod_phi(P)> child) : x(value), p(child) {
    }

    // 現在の法の正規化された剰余を返す。計算量 O(1)。
    constexpr Int val() const {
        return x < P ? x : x - P;
    }

    // φ鎖の各法で和を求める。計算量 O(L)。
    friend constexpr expmodint operator+(expmodint a, expmodint b) {
        return {detail::expmod_getmod<P>(a.x + b.x), a.p + b.p};
    }

    // φ鎖の各法で積を求める。計算量 O(L)。
    friend constexpr expmodint operator*(expmodint a, expmodint b) {
        return {detail::expmod_mul<P>(a.x, b.x), a.p * b.p};
    }

    // 現在の法での値をストリームへ出力する。計算量 O(log P)。
    friend std::ostream &operator<<(std::ostream &s, expmodint a) {
        return s << a.val();
    }
};

template <> struct expmodint<1> {
    Int x = 0;
    // 法1の零を構築する。計算量 O(1)。
    constexpr expmodint() = default;

    // 法1の剰余を構築する。計算量 O(1)。
    constexpr explicit expmodint(Int) {
    }

    // 法1の剰余0を返す。計算量 O(1)。
    constexpr Int val() const {
        return 0;
    }

    // 法1で加算する。計算量 O(1)。
    friend constexpr expmodint operator+(expmodint, expmodint) {
        return {};
    }

    // 法1で乗算する。計算量 O(1)。
    friend constexpr expmodint operator*(expmodint, expmodint) {
        return {};
    }

    // 法1の剰余を出力する。計算量 O(1)。
    friend std::ostream &operator<<(std::ostream &s, expmodint) {
        return s << 0;
    }
};

// 整数から expmodint を作る。φ鎖の長さを L として O(L)。
template <Int P> constexpr auto tomodint(Int x) {
    return expmodint<P>(x);
}

// 整数から expmodint を作る。φ鎖の長さを L として O(L)。
template <Int P> constexpr auto toModint(Int x) {
    return expmodint<P>(x);
}

// 現在の法での剰余を返す。計算量 O(1)。
template <Int P> constexpr Int val(expmodint<P> a) {
    return a.val();
}

// φ鎖を用いて累乗を求める。鎖の長さを L、指数の表現値の上限を B として O(L log(B+1))。
template <Int P, Int Q> constexpr expmodint<P> pow(expmodint<P> a, expmodint<Q> b) {
    if constexpr (P == 1)
        return {};
    else
        return {detail::expmod_pow<P>(a.x, b.p.x), pow(a.p, b.p)};
}

// 剰余を十進文字列に変換する。計算量 O(log(P+1))。
template <Int P> std::string to_string(expmodint<P> a) {
    return std::to_string(a.val());
}
}
