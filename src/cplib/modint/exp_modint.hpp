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

// 元実装のphi鎖と[P,2P)の印を保持する版。積がIntを越えない入力を前提とする。
template <Int P> struct expmodint {
    static_assert(P > 1);
    Int x = 0;
    expmodint<detail::expmod_phi(P)> p;
    constexpr expmodint() = default;

    constexpr explicit expmodint(Int value) : x(detail::expmod_getmod<P>(value)), p(value) {
    }

    constexpr expmodint(Int value, expmodint<detail::expmod_phi(P)> child) : x(value), p(child) {
    }

    constexpr Int val() const {
        return x < P ? x : x - P;
    }

    friend constexpr expmodint operator+(expmodint a, expmodint b) {
        return {detail::expmod_getmod<P>(a.x + b.x), a.p + b.p};
    }

    friend constexpr expmodint operator*(expmodint a, expmodint b) {
        return {detail::expmod_mul<P>(a.x, b.x), a.p * b.p};
    }

    friend std::ostream &operator<<(std::ostream &s, expmodint a) {
        return s << a.val();
    }
};

template <> struct expmodint<1> {
    Int x = 0;
    constexpr expmodint() = default;

    constexpr explicit expmodint(Int) {
    }

    constexpr Int val() const {
        return 0;
    }

    friend constexpr expmodint operator+(expmodint, expmodint) {
        return {};
    }

    friend constexpr expmodint operator*(expmodint, expmodint) {
        return {};
    }

    friend std::ostream &operator<<(std::ostream &s, expmodint) {
        return s << 0;
    }
};

template <Int P> constexpr auto tomodint(Int x) {
    return expmodint<P>(x);
}

template <Int P> constexpr auto toModint(Int x) {
    return expmodint<P>(x);
}

template <Int P> constexpr Int val(expmodint<P> a) {
    return a.val();
}

template <Int P, Int Q> constexpr expmodint<P> pow(expmodint<P> a, expmodint<Q> b) {
    if constexpr (P == 1)
        return {};
    else
        return {detail::expmod_pow<P>(a.x, b.p.x), pow(a.p, b.p)};
}

template <Int P> std::string to_string(expmodint<P> a) {
    return std::to_string(a.val());
}
}
