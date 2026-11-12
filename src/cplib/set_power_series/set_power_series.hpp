#pragma once
#include <cplib/convolution/subset_convolution.hpp>
#include <cplib/modint/modint.hpp>
#include <cplib/math/isprime.hpp>

namespace cplib {
namespace detail {
inline Int setVariableCount(std::size_t size) {
    assert(size > 0 && std::has_single_bit(size));
    return std::countr_zero(size);
}
template <class T>
std::vector<T> setComposeDerivatives(const std::vector<T> &derivatives, const std::vector<T> &s) {
    Int n = setVariableCount(s.size());
    std::vector<std::vector<T>> values(n + 1);
    for (Int k = 0; k <= n; ++k)
        values[k] = {derivatives[k]};
    std::size_t width = 1;
    for (Int v = 0; v < n; ++v) {
        std::vector<T> h(s.begin() + width, s.begin() + width * 2);
        for (Int k = 0; k < n - v; ++k) {
            auto upper = subsetConvolution(values[k + 1], h);
            values[k].insert(values[k].end(), upper.begin(), upper.end());
        }
        values.resize(n - v);
        width *= 2;
    }
    return std::move(values[0]);
}
template <class T> T setUnitInverse(T c) {
    assert(c != T(0));
    T inverse = c.inv();
    assert(c * inverse == T(1));
    return inverse;
}
template <class T> T setScalarPower(T c, UInt exponent) {
    T result = 1;
    while (exponent > 0) {
        if (exponent & 1)
            result *= c;
        exponent >>= 1;
        if (exponent > 0)
            c *= c;
    }
    return result;
}
}
// exp(s)を返す。長さ2^n、s[0]=0、素数法>nが必要。時間 O(n²2^n)、領域 O(n2^n)。
template <class T> std::vector<T> setExp(const std::vector<T> &s) {
    Int n = detail::setVariableCount(s.size());
    assert(s[0] == T(0) && T::umod() > n && isprime(Int(T::umod())));
    std::vector<T> result{T(1)};
    for (std::size_t width = 1; width < s.size(); width *= 2) {
        auto upper = subsetConvolution(result, std::vector<T>(s.begin() + width, s.begin() + width * 2));
        result.insert(result.end(), upper.begin(), upper.end());
    }
    return result;
}
// log(s)を返す。s[0]=1、素数法>nが必要。時間 O(n²2^n)、領域 O(n2^n)。
template <class T> std::vector<T> setLog(const std::vector<T> &s) {
    Int n = detail::setVariableCount(s.size());
    assert(s[0] == T(1) && T::umod() > n && isprime(Int(T::umod())));
    std::vector<T> derivatives(n + 1);
    if (n > 0)
        derivatives[1] = 1;
    for (Int k = 2; k <= n; ++k)
        derivatives[k] = derivatives[k - 1] * (-(k - 1));
    return detail::setComposeDerivatives(derivatives, s);
}
// 定数項が単元のsの逆元を返す。時間 O(n²2^n)、領域 O(n2^n)。
template <class T> std::vector<T> setInv(const std::vector<T> &s) {
    Int n = detail::setVariableCount(s.size());
    T inverse = detail::setUnitInverse(s[0]);
    std::vector<T> derivatives(n + 1);
    derivatives[0] = inverse;
    for (Int k = 1; k <= n; ++k)
        derivatives[k] = derivatives[k - 1] * (-k) * inverse;
    return detail::setComposeDerivatives(derivatives, s);
}
// sの整数冪を返す。負指数では定数項が単元であること。時間 O(n²2^n+n log(|e|+1))、領域 O(n2^n)。
template <class T> std::vector<T> setPow(const std::vector<T> &s, Int exponent) {
    Int n = detail::setVariableCount(s.size());
    std::vector<T> derivatives(n + 1);
    if (exponent >= 0) {
        T falling = 1;
        for (Int k = 0; k <= std::min(n, exponent); ++k) {
            derivatives[k] = falling * detail::setScalarPower(s[0], UInt(exponent - k));
            falling *= T(exponent) - T(k);
        }
    } else {
        T inverse = detail::setUnitInverse(s[0]);
        derivatives[0] = detail::setScalarPower(inverse, UInt(-(exponent + 1)) + 1);
        for (Int k = 1; k <= n; ++k)
            derivatives[k] = derivatives[k - 1] * (T(exponent) - T(k - 1)) * inverse;
    }
    return detail::setComposeDerivatives(derivatives, s);
}
// 定数項の平方根を指定して平方根を返す。2と定数項は単元が必要。全零はroot=0も可。
// 時間 O(n²2^n)、領域 O(n2^n)。
template <class T> std::vector<T> setSqrt(const std::vector<T> &s, T constantRoot) {
    Int n = detail::setVariableCount(s.size());
    assert(constantRoot * constantRoot == s[0]);
    if (constantRoot == T(0)) {
        for (T c : s)
            assert(c == T(0));
        return std::vector<T>(s.size());
    }
    T inverse = detail::setUnitInverse(s[0]), half = detail::setUnitInverse(T(2));
    std::vector<T> derivatives(n + 1);
    derivatives[0] = constantRoot;
    for (Int k = 1; k <= n; ++k)
        derivatives[k] = derivatives[k - 1] * (half - T(k - 1)) * inverse;
    return detail::setComposeDerivatives(derivatives, s);
}
// 多項式f(s)を返す。M=f.size()として時間 O(n²2^n+nM)、領域 O(n2^n+M)。
template <class T>
std::vector<T> polynomialCompositeSetPowerSeries(std::vector<T> f, const std::vector<T> &s) {
    Int n = detail::setVariableCount(s.size());
    std::vector<T> derivatives(n + 1);
    for (Int k = 0; k <= n && !f.empty(); ++k) {
        for (std::size_t j = f.size(); j-- > 0;)
            derivatives[k] = derivatives[k] * s[0] + f[j];
        for (std::size_t j = 1; j < f.size(); ++j)
            f[j - 1] = f[j] * T(j);
        f.pop_back();
    }
    return detail::setComposeDerivatives(derivatives, s);
}
// Σ_S weights[S](s^k)[S]をk=0..m-1について返す。時間 O(n²2^n+nm)、領域 O(n2^n+m)。
template <class T>
std::vector<T> setPowerProjection(const std::vector<T> &s, const std::vector<T> &weights, Int m) {
    Int n = detail::setVariableCount(s.size());
    assert(s.size() == weights.size() && m >= 0);
    if (m == 0)
        return {};
    std::vector<std::vector<T>> adjoints{weights};
    for (Int v = n; v-- > 0;) {
        std::size_t width = std::size_t(1) << v;
        std::vector<T> h(s.begin() + width, s.begin() + width * 2);
        std::vector<std::vector<T>> next(adjoints.size() + 1, std::vector<T>(width));
        for (std::size_t k = 0; k < adjoints.size(); ++k) {
            std::vector<T> reversed(width);
            for (std::size_t mask = 0; mask < width; ++mask) {
                next[k][mask] += adjoints[k][mask];
                reversed[width - 1 - mask] = adjoints[k][width + mask];
            }
            auto product = subsetConvolution(h, reversed);
            for (std::size_t mask = 0; mask < width; ++mask)
                next[k + 1][mask] += product[width - 1 - mask];
        }
        adjoints = std::move(next);
    }
    std::vector<T> powers(m), result(m);
    powers[0] = 1;
    for (Int k = 1; k < m; ++k)
        powers[k] = powers[k - 1] * s[0];
    for (Int k = 0; k < m; ++k) {
        T falling = 1;
        for (Int j = 0; j <= std::min(k, n); ++j) {
            result[k] += falling * powers[k - j] * adjoints[j][0];
            falling *= T(k) - T(j);
        }
    }
    return result;
}
}
