#pragma once
#include <cplib/common.hpp>
#include <cplib/math/inv_gcd.hpp>

namespace cplib::detail {
template <class T>
std::vector<T> convolution_schoolbook(const std::vector<T> &f, const std::vector<T> &g) {
    if (f.empty() || g.empty())
        return {};
    std::vector<T> out(f.size() + g.size() - 1, T(0));
    if (f.size() > g.size()) {
        for (std::size_t i = 0; i < f.size(); ++i)
            for (std::size_t j = 0; j < g.size(); ++j)
                out[i + j] += f[i] * g[j];
    } else {
        for (std::size_t j = 0; j < g.size(); ++j)
            for (std::size_t i = 0; i < f.size(); ++i)
                out[i + j] += f[i] * g[j];
    }
    return out;
}

inline Int convolution_floor_mod(Int x, Int m) {
    Int r = x % m;
    return r < 0 ? r + m : r;
}

inline std::vector<Int> convolution_crt_ll(const std::vector<std::uint32_t> &c1,
                                           const std::vector<std::uint32_t> &c2,
                                           const std::vector<std::uint32_t> &c3) {
    constexpr UInt M1 = 754974721, M2 = 167772161, M3 = 469762049, M12 = M1 * M2, M23 = M2 * M3,
                   M31 = M3 * M1, M123 = M1 * M2 * M3;
    static const UInt i1 = inv_gcd(M23, M1).second, i2 = inv_gcd(M31, M2).second,
                      i3 = inv_gcd(M12, M3).second;
    constexpr UInt offset[] = {0, 0, M123, 2 * M123, 3 * M123};
    std::vector<Int> out(c1.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        UInt x = UInt(c1[i]) * i1 % M1 * M23;
        x += UInt(c2[i]) * i2 % M2 * M31;
        x += UInt(c3[i]) * i3 % M3 * M12;
        Int diff = Int(c1[i]) - convolution_floor_mod(std::bit_cast<Int>(x), M1);
        if (diff < 0)
            diff += M1;
        x -= offset[diff % 5];
        out[i] = std::bit_cast<Int>(x);
    }
    return out;
}
}
