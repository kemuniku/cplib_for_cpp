#pragma once
#include <cplib/fps/formal_power_series.hpp>

namespace cplib {
namespace detail {
template <class T> std::vector<T> parityTerms(const std::vector<T> &f, Int parity) {
    if (Int(f.size()) <= parity)
        return {};
    std::vector<T> out;
    out.reserve((f.size() - parity + 1) / 2);
    for (Int i = parity; i < Int(f.size()); i += 2)
        out.push_back(f[i]);
    return out;
}
}

// 有理型母関数の第k係数。O(M(D)log(k+1))。
template <Modint T>
T bostanMori(const std::vector<T> &numerator, const std::vector<T> &denominator, Int k) {
    assert(k >= 0 && !denominator.empty() && denominator[0].val() != 0);
    auto p = normalized(numerator), q = normalized(denominator);
    while (k > 0) {
        auto neg = q;
        for (std::size_t i = 1; i < neg.size(); i += 2)
            neg[i] = -neg[i];
        p = detail::parityTerms(p * neg, k & 1);
        q = detail::parityTerms(q * neg, 0);
        k >>= 1;
    }
    return p.empty() ? T(0) : p[0] / q[0];
}

// a[n] = sum(coefficients[i] * a[n-i-1], i=0..<d) の第k項を求める。
template <Modint T>
T linearRecurrenceKth(const std::vector<T> &initial, const std::vector<T> &coefficients, Int k) {
    assert(initial.size() == coefficients.size() && !initial.empty() && k >= 0);
    if (k < Int(initial.size()))
        return initial[k];
    std::vector<T> q(coefficients.size() + 1);
    q[0] = 1;
    for (std::size_t i = 0; i < coefficients.size(); ++i)
        q[i + 1] = -coefficients[i];
    return bostanMori(prefix(initial * q, coefficients.size()), q, k);
}
}
