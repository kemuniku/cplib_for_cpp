#pragma once
#include <cplib/modint/barrett_impl.hpp>
#include <cplib/modint/montgomery_impl.hpp>
#include <cplib/math/isqrt.hpp>
#include <tuple>

namespace cplib {
using modint998244353_montgomery = StaticMontgomeryModint<998244353>;
using modint1000000007_montgomery = StaticMontgomeryModint<1000000007>;
using modint_montgomery = DynamicMontgomeryModint<1>;
using modint998244353_barrett = StaticBarrettModint<998244353>;
using modint1000000007_barrett = StaticBarrettModint<1000000007>;
using modint_barrett = DynamicBarrettModint<1>;

// 分母1..ubの既約分数を列挙し、|分子|+分母が最小のものを返す。O(ub log ub)。
template <Modint M> std::string estimate_rational(M a, Int ub) {
    assert(ub >= 1);
    std::vector<std::tuple<Int, Int, Int>> v;
    v.reserve(ub);
    for (Int d = 1; d <= ub; ++d) {
        Int n = (a * d).val();
        if (n * 2 > M::mod())
            n = -(M::mod() - n);
        if (std::gcd(n, d) > 1)
            continue;
        v.emplace_back((n < 0 ? -n : n) + d, n, d);
    }
    std::sort(v.begin(), v.end());
    return std::to_string(std::get<1>(v[0])) + "/" + std::to_string(std::get<2>(v[0]));
}

template <Modint M> std::string estimate_rational(M a) {
    return estimate_rational(a, isqrt(Int(M::mod())));
}

// Nimのstatic int除算に相当。静的modintなら逆元を定数評価する。
template <Int Divisor, Modint M> constexpr M divide_constant(M a) {
    if constexpr (requires { std::integral_constant<Int, M::mod()>{}; }) {
        constexpr M inverse = M(Divisor).inv();
        return a * inverse;
    } else
        return a * M(Divisor).inv();
}
}
