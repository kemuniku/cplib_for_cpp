#pragma once
#include <cplib/fps/formal_power_series.hpp>

namespace cplib {
namespace detail {
// Kinoshita–Li法。二変数積を平坦化し、転置した冪射影を再帰的に復元する。
template <Modint T>
std::vector<std::vector<T>> compositionRec(const std::vector<T> &outer,
                                           const std::vector<std::vector<T>> &denominator, Int n,
                                           Int yDegree) {
    if (n == 0) {
        std::vector<std::vector<T>> output_value(yDegree, std::vector<T>(1));
        for (Int i = 0; i < std::min(Int(outer.size()), yDegree); ++i)
            output_value[yDegree - 1 - i][0] = outer[i];
        return output_value;
    }
    Int degree = 0;
    for (const auto &p : denominator)
        degree = std::max(degree, std::min(Int(p.size()), n + 1) - 1);
    Int stride = n + degree + 1, length = yDegree * stride + degree + 1;
    std::vector<T> positive(length), negative(length);
    for (Int y = 0; y <= yDegree; ++y)
        for (Int x = 0; x < std::min(Int(denominator[y].size()), n + 1); ++x) {
            positive[y * stride + x] = denominator[y][x];
            negative[y * stride + x] = (x & 1) ? -denominator[y][x] : denominator[y][x];
        }
    auto product = positive * negative;
    Int half = n / 2;
    std::vector<std::vector<T>> next(2 * yDegree + 1, std::vector<T>(half + 1));
    for (Int y = 0; y <= 2 * yDegree; ++y)
        for (Int x = 0; x <= half; ++x) {
            Int i = y * stride + 2 * x;
            if (i < Int(product.size()))
                next[y][x] = product[i];
        }
    auto projected = compositionRec(outer, next, half, 2 * yDegree);
    std::vector<T> lifted((2 * yDegree - 1) * stride + 2 * half + 1);
    for (Int y = 0; y < 2 * yDegree; ++y)
        for (Int x = 0; x < Int(projected[y].size()); ++x)
            lifted[y * stride + 2 * x] = projected[y][x];
    auto recovered = lifted * negative;
    std::vector<std::vector<T>> output_value(yDegree, std::vector<T>(n + 1));
    for (Int y = 0; y < yDegree; ++y)
        for (Int x = 0; x <= n; ++x) {
            Int i = (yDegree + y) * stride + x;
            if (i < Int(recovered.size()))
                output_value[y][x] = recovered[i];
        }
    return output_value;
}
}

// outer(inner(x)) mod x^n。inner(0)=0。O(n log²n)。
template <Modint T>
std::vector<T> compose(const std::vector<T> &outer, const std::vector<T> &inner, Int n) {
    if (n <= 0)
        return {};
    assert(inner.empty() || inner[0].val() == 0);
    if (n == 1)
        return {outer.empty() ? T(0) : outer[0]};
    std::vector<std::vector<T>> denominator = {{T(1)}, -prefix(inner, n)};
    return prefix(detail::compositionRec(prefix(outer, n), denominator, n - 1, 1)[0], n);
}

// outer.size() 項まで合成する。計算量は compose(outer, inner, outer.size()) と同じ。
template <Modint T> auto compose(const std::vector<T> &outer, const std::vector<T> &inner) {
    return compose(outer, inner, outer.size());
}

// Newton法。高次誤差のみを導関数比で補正する。
// Newton法により f(g(x)) = x (mod x^n) を満たすgを求める。NTT使用時 O(n log²(n+1))。
template <Modint T> std::vector<T> compositionalInverse(const std::vector<T> &f, Int n) {
    if (n <= 0)
        return {};
    assert(f.size() >= 2 && f[0].val() == 0 && f[1].val() != 0);
    if (n == 1)
        return std::vector<T>(1);
    std::vector<T> output_value = {T(0), f[1].inv()};
    for (Int m = 2; m < n;) {
        Int next = std::min(m * 2, n);
        auto composed = compose(f, output_value, next);
        Int size = next - m;
        std::vector<T> error(composed.begin() + m, composed.begin() + next);
        auto inverseSlope = prefix(derivative(output_value) * inv(derivative(composed), size), size);
        auto correction = prefix(error * inverseSlope, size);
        output_value.resize(next);
        for (Int i = 0; i < size; ++i)
            output_value[m + i] -= correction[i];
        m = next;
    }
    return output_value;
}

// f.size() 項まで合成逆元を求める。計算量は compositionalInverse(f, f.size()) と同じ。
template <Modint T> auto compositionalInverse(const std::vector<T> &f) {
    return compositionalInverse(f, f.size());
}
}
