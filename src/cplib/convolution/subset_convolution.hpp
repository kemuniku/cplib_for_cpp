#pragma once
#include <cplib/convolution/set_transform.hpp>
#include <cplib/convolution/private/subset_convolution_barrett.hpp>

namespace cplib {
namespace detail {
template <class T> struct SubsetBarrettModulus {
    static constexpr std::uint32_t value = 0;
};
template <std::uint32_t M> struct SubsetBarrettModulus<BasicBarrettModint<M, false>> {
    static constexpr std::uint32_t value = M;
};
}

// c[S]=Σ(U⊆S) a[U]b[S xor U]を返す。同じ2冪長、両方空なら空。
// 長さ2^nとして時間 O(n² 2^n)、領域 O(n 2^n)。
template <class T>
std::vector<T> subsetConvolution(const std::vector<T> &a, const std::vector<T> &b) {
    assert(a.size() == b.size());
    if (a.empty())
        return {};
    assert(std::has_single_bit(a.size()));
    const std::size_t size = a.size();
    const Int n = std::countr_zero(size);
    constexpr auto modulus = detail::SubsetBarrettModulus<T>::value;
    if constexpr (modulus > 1 && modulus <= (1u << 30)) {
        if (n <= 20)
            return subsetConvolutionBarrett(a, b);
    }
    const Int stride = n + 1;
    std::vector<T> left(size * stride), right(size * stride);
    for (std::size_t mask = 0; mask < size; ++mask) {
        Int rank = std::popcount(mask);
        left[mask * stride + rank] = a[mask];
        right[mask * stride + rank] = b[mask];
    }
    for (std::size_t bit = 1; bit < size; bit *= 2) {
        for (std::size_t mask = 0; mask < size; ++mask) {
            if (!(mask & bit))
                continue;
            for (Int k = 0; k <= n; ++k) {
                left[mask * stride + k] += left[(mask ^ bit) * stride + k];
                right[mask * stride + k] += right[(mask ^ bit) * stride + k];
            }
        }
    }
    for (std::size_t mask = 0; mask < size; ++mask) {
        const auto offset = mask * stride;
        for (Int k = n; k >= 0; --k) {
            T sum{};
            for (Int j = 0; j <= k; ++j)
                sum += left[offset + j] * right[offset + k - j];
            left[offset + k] = sum;
        }
    }
    for (std::size_t bit = 1; bit < size; bit *= 2)
        for (std::size_t mask = 0; mask < size; ++mask)
            if (mask & bit)
                for (Int k = 0; k <= n; ++k)
                    left[mask * stride + k] -= left[(mask ^ bit) * stride + k];
    std::vector<T> result(size);
    for (std::size_t mask = 0; mask < size; ++mask)
        result[mask] = left[mask * stride + std::popcount(mask)];
    return result;
}
}
