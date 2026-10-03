#pragma once
#include <cplib/modint/modint.hpp>
#include <cplib/math/isprime.hpp>
#include <cplib/convolution/private/convolution_common.hpp>
#include <cplib/convolution/private/avx2_ntt.hpp>

namespace cplib {
namespace detail {
inline std::pair<std::uint32_t, bool> nttPrimalityCache{};

inline bool isNttFriendlyModulus(std::uint32_t modulus, std::uint32_t size) {
    if (modulus <= 1 || modulus >= (1u << 30) || (modulus - 1) % size)
        return false;
    if (nttPrimalityCache.first != modulus)
        nttPrimalityCache = {modulus, isprime(Int(modulus))};
    return nttPrimalityCache.second;
}

inline std::vector<std::uint32_t> convolutionNttFriendlyU32(const std::vector<std::uint32_t> &f,
                                                            const std::vector<std::uint32_t> &g,
                                                            std::uint32_t modulus,
                                                            std::uint32_t primitiveRoot) {
    if (f.empty() || g.empty())
        return {};
    std::size_t deg = f.size() + g.size() - 1;
    if (std::min(f.size(), g.size()) <= 60) {
        std::vector<std::uint32_t> out(deg);
        if (f.size() > g.size()) {
            for (std::size_t i = 0; i < f.size(); ++i)
                for (std::size_t j = 0; j < g.size(); ++j)
                    out[i + j] = (out[i + j] + UInt(f[i]) * g[j]) % modulus;
        } else {
            for (std::size_t j = 0; j < g.size(); ++j)
                for (std::size_t i = 0; i < f.size(); ++i)
                    out[i + j] = (out[i + j] + UInt(f[i]) * g[j]) % modulus;
        }
        return out;
    }
    std::size_t size = std::bit_ceil(deg);
    std::vector<std::uint32_t> out(size);
    avx2_ntt::convolution_ntt_friendly(out.data(), f.data(), f.size(), g.data(), g.size(), size,
                                       modulus, primitiveRoot, false);
    out.resize(deg);
    return out;
}

template <class T>
std::vector<T> convolutionArbitraryMod(const std::vector<T> &f, const std::vector<T> &g) {
    constexpr UInt M1 = 754974721, M2 = 167772161, M3 = 469762049, M12 = M1 * M2;
    static const UInt inv1 = inv_gcd(M1 % M2, M2).second, inv12 = inv_gcd(M12 % M3, M3).second;
    UInt target = T::umod();
    assert(target > 0 && target < (UInt(1) << 31));
    assert(std::bit_ceil(f.size() + g.size() - 1) <= (std::size_t(1) << 24));
    std::vector<std::uint32_t> fm(f.size()), gm(g.size());
    auto run = [&](UInt modulus, std::uint32_t root) {
        for (std::size_t i = 0; i < f.size(); ++i)
            fm[i] = UInt(f[i].val()) % modulus;
        for (std::size_t i = 0; i < g.size(); ++i)
            gm[i] = UInt(g[i].val()) % modulus;
        return convolutionNttFriendlyU32(fm, gm, modulus, root);
    };
    auto c1 = run(M1, 11), c2 = run(M2, 3), c3 = run(M3, 3);
    UInt m1Target = M1 % target, m12Target = M12 % target;
    std::vector<T> out(c1.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        UInt r1 = c1[i], t2 = (UInt(c2[i]) + M2 - r1 % M2) % M2 * inv1 % M2,
             r12 = (r1 + (M1 % M3) * t2) % M3, t3 = (UInt(c3[i]) + M3 - r12) % M3 * inv12 % M3;
        out[i] = T((r1 % target + m1Target * t2 % target + m12Target * t3 % target) % target);
    }
    return out;
}

// memcpy経由でオブジェクト表現を受け渡し、strict aliasing違反を避ける。
template <class T> auto modintRawWords(const std::vector<T> &f) {
    static_assert(sizeof(T) == sizeof(std::uint32_t) && std::is_trivially_copyable_v<T>);
    std::vector<std::uint32_t> out(f.size());
    if (!f.empty())
        std::memcpy(out.data(), f.data(), sizeof(T) * f.size());
    return out;
}

template <class T> auto modintFromRawWords(const std::vector<std::uint32_t> &f, std::size_t size) {
    std::vector<T> out(size);
    for (std::size_t i = 0; i < size; ++i)
        out[i] = std::bit_cast<T>(f[i]);
    return out;
}
}

template <class T> auto convolution_naive(const std::vector<T> &f, const std::vector<T> &g) {
    return detail::convolution_schoolbook(f, g);
}

// 現行版：AVX2の半分ゼロ特化・radix-4 NTT、任意法は3素数CRT。O(N log N)。
template <class T>
    requires(BarrettModint<T> || MontgomeryModint<T>)
std::vector<T> convolution(const std::vector<T> &f, const std::vector<T> &g) {
    if (f.empty() || g.empty())
        return {};
    if (std::min(f.size(), g.size()) <= 60)
        return convolution_naive(f, g);
    std::size_t deg = f.size() + g.size() - 1, size = std::bit_ceil(deg);
    if (detail::isNttFriendlyModulus(T::umod(), size)) {
        auto a = detail::modintRawWords(f), b = detail::modintRawWords(g);
        std::vector<std::uint32_t> out(size);
        detail::avx2_ntt::convolution_ntt_friendly(out.data(), a.data(), a.size(), b.data(),
                                                   b.size(), size, T::umod(), 0,
                                                   MontgomeryModint<T>);
        return detail::modintFromRawWords<T>(out, deg);
    }
    return detail::convolutionArbitraryMod(f, g);
}

// 2冪長の巡回積。Montgomery版は通常表現へ変換してカーネルを呼ぶ元経路を保持。
template <class T>
    requires(BarrettModint<T> || MontgomeryModint<T>)
std::vector<T> convolutionCyclicPowerOfTwo(const std::vector<T> &f, const std::vector<T> &g,
                                           Int n) {
    assert(n > 0 && std::has_single_bit(UInt(n)) && f.size() <= std::size_t(n) &&
           g.size() <= std::size_t(n));
    std::vector<T> out(n);
    if (f.empty() || g.empty())
        return out;
    if (n >= 64 && detail::isNttFriendlyModulus(T::umod(), n)) {
        std::vector<std::uint32_t> a(f.size()), b(g.size()), res(n);
        if constexpr (MontgomeryModint<T>) {
            for (std::size_t i = 0; i < f.size(); ++i)
                a[i] = f[i].val();
            for (std::size_t i = 0; i < g.size(); ++i)
                b[i] = g[i].val();
        } else {
            a = detail::modintRawWords(f);
            b = detail::modintRawWords(g);
        }
        detail::avx2_ntt::convolution_ntt_friendly(res.data(), a.data(), a.size(), b.data(),
                                                   b.size(), n, T::umod(), 0, false);
        if constexpr (MontgomeryModint<T>) {
            for (Int i = 0; i < n; ++i)
                out[i] = T(res[i]);
        } else
            out = detail::modintFromRawWords<T>(res, n);
        return out;
    }
    auto product = convolution(f, g);
    for (std::size_t i = 0; i < product.size(); ++i)
        out[i < std::size_t(n) ? i : i - n] += product[i];
    return out;
}

template <Int M>
std::vector<Int> convolution(const std::vector<Int> &f, const std::vector<Int> &g) {
    static_assert(M > 0 && M < (Int(1) << 31));
    using Mint = StaticBarrettModint<std::uint32_t(M)>;
    std::vector<Mint> a(f.begin(), f.end()), b(g.begin(), g.end());
    auto product = convolution(a, b);
    std::vector<Int> out(product.size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = product[i].val();
    return out;
}

// 各係数が符号付き64ビットに収まる整数畳み込み。CRTの2^64桁あふれ補正を保持。
inline std::vector<Int> convolution_ll(const std::vector<Int> &f, const std::vector<Int> &g) {
    if (f.empty() || g.empty())
        return {};
    std::vector<std::uint32_t> a(f.size()), b(g.size());
    auto run = [&](Int modulus, std::uint32_t root) {
        for (std::size_t i = 0; i < f.size(); ++i)
            a[i] = detail::convolution_floor_mod(f[i], modulus);
        for (std::size_t i = 0; i < g.size(); ++i)
            b[i] = detail::convolution_floor_mod(g[i], modulus);
        return detail::convolutionNttFriendlyU32(a, b, modulus, root);
    };
    auto c1 = run(754974721, 11), c2 = run(167772161, 3), c3 = run(469762049, 3);
    return detail::convolution_crt_ll(c1, c2, c3);
}
}
