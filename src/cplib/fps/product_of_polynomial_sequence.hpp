#pragma once
#include <cplib/convolution/convolution.hpp>

namespace cplib {
// 次数の和で分割位置を均衡させる多項式列積。Barrett/998244353は専用AVX2カーネル。
template <Modint T>
std::vector<T> productOfPolynomialSequence(const std::vector<std::vector<T>> &polynomials) {
    if (polynomials.empty())
        return {T(1)};
    Int totalLength = 1;
    T scalar = 1;
    std::vector<const std::vector<T> *> factors;
    factors.reserve(polynomials.size());
    for (const auto &p : polynomials) {
        if (p.empty())
            return {};
        totalLength += p.size() - 1;
        if (p.size() == 1)
            scalar *= p[0];
        else
            factors.push_back(&p);
    }
    if (scalar.val() == 0)
        return std::vector<T>(totalLength);
    if (factors.empty())
        return {scalar};
    if constexpr (BarrettModint<T>)
        if (T::umod() == 998244353) {
            std::vector<std::vector<std::uint32_t>> words;
            words.reserve(factors.size());
            std::vector<const std::uint32_t *> pointers;
            std::vector<std::size_t> sizes;
            for (auto p : factors) {
                words.push_back(detail::modintRawWords(*p));
                pointers.push_back(words.back().data());
                sizes.push_back(p->size());
            }
            std::vector<std::uint32_t> output(totalLength);
            detail::avx2_ntt::product_polynomial_sequence_998(output.data(), pointers.data(),
                                                              sizes.data(), factors.size());
            auto out = detail::modintFromRawWords<T>(output, totalLength);
            if (scalar.val() != 1)
                for (auto &x : out)
                    x *= scalar;
            return out;
        }
    std::vector<T> out;
    if (factors.size() == 1)
        out = *factors[0];
    else {
        std::vector<Int> degrees(factors.size() + 1);
        for (Int i = 0; i < Int(factors.size()); ++i)
            degrees[i + 1] = degrees[i] + factors[i]->size() - 1;
        auto middle = [&](Int left, Int right) {
            if (right - left == 2)
                return left + 1;
            Int target = degrees[left] + (degrees[right] - degrees[left]) / 2, low = left + 1,
                high = right;
            while (low < high) {
                Int m = (low + high) / 2;
                if (degrees[m] < target)
                    low = m + 1;
                else
                    high = m;
            }
            Int result = low;
            if (result > left + 1) {
                Int current = std::abs((degrees[result] - degrees[left]) -
                                       (degrees[right] - degrees[result])),
                    previous = std::abs((degrees[result - 1] - degrees[left]) -
                                        (degrees[right] - degrees[result - 1]));
                if (previous < current)
                    --result;
            }
            return result;
        };
        auto solve = [&](auto &&self, Int left, Int right) -> std::vector<T> {
            if (left + 1 == right)
                return *factors[left];
            Int m = middle(left, right);
            auto a = self(self, left, m), b = self(self, m, right);
            return convolution(a, b);
        };
        out = solve(solve, 0, factors.size());
    }
    if (scalar.val() != 1)
        for (auto &x : out)
            x *= scalar;
    return out;
}
}
