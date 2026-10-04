#pragma once
#include <cplib/common.hpp>
#include <type_traits>

namespace cplib {
// 正の2冪長のアダマール変換。O(N log N)。
template <class T> void FastHadamardTransForm(std::vector<T> &u) {
    std::size_t n = u.size();
    assert(n > 0 && std::has_single_bit(n));
    for (std::size_t i = 1; i < n; i <<= 1)
        for (std::size_t j = 0; j < n; ++j)
            if (!(j & i)) {
                T x = u[j], y = u[j + i];
                u[j] = x + y;
                u[j + i] = x - y;
            }
}

// 同じ正の2冪長の配列のXOR畳み込みを返す。O(N log N)。
template <class T> std::vector<T> xorConvolution(std::vector<T> u, std::vector<T> v) {
    assert(u.size() == v.size());
    FastHadamardTransForm(u);
    FastHadamardTransForm(v);
    for (std::size_t i = 0; i < u.size(); ++i)
        u[i] *= v[i];
    FastHadamardTransForm(u);
    if constexpr (std::is_same_v<T, Int>) {
        int k = std::countr_zero(u.size());
        for (auto &x : u)
            x >>= k;
    } else {
        T inv = T(1) / T(u.size());
        for (auto &x : u)
            x *= inv;
    }
    return u;
}
}
