#pragma once
// CPLIB_HAMMING_DISTANCE_SCALARで通常版、CPLIB_HAMMING_DISTANCE_AVX2でAVX-512を無効化する。
#include <cplib/common.hpp>
#include <string_view>
#include <stdexcept>
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
namespace cplib::detail {
__attribute__((target("avx2,popcnt")))
inline std::size_t cplib_hamming_avx2(const unsigned char *s, const unsigned char *t, size_t n) {
    // 32バイトずつ比較して不一致数を加算する。O(N)。
    size_t i = 0, count = 0;
    for (; n - i >= 32; i += 32) {
        __m256i a = _mm256_loadu_si256((const __m256i *)(s + i));
        __m256i b = _mm256_loadu_si256((const __m256i *)(t + i));
        unsigned int mask = (unsigned int)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b));
        count += 32 - __builtin_popcount(mask);
    }
    for (; i < n; ++i) count += s[i] != t[i];
    return count;
}

__attribute__((target("avx512f,avx512bw,popcnt")))
inline std::size_t cplib_hamming_avx512(const unsigned char *s, const unsigned char *t, size_t n) {
    // 64バイトずつ比較し、端数はマスク付きで読む。O(N)。
    size_t i = 0, count = 0;
    for (; n - i >= 64; i += 64) {
        __m512i a = _mm512_loadu_si512((const void *)(s + i));
        __m512i b = _mm512_loadu_si512((const void *)(t + i));
        count += __builtin_popcountll((unsigned long long)_mm512_cmpneq_epi8_mask(a, b));
    }
    if (i < n) {
        __mmask64 valid = (__mmask64)((UINT64_C(1) << (n - i)) - 1);
        __m512i a = _mm512_maskz_loadu_epi8(valid, (const void *)(s + i));
        __m512i b = _mm512_maskz_loadu_epi8(valid, (const void *)(t + i));
        count += __builtin_popcountll((unsigned long long)(_mm512_cmpneq_epi8_mask(a, b) & valid));
    }
    return count;
}

}
#endif
namespace cplib {
// バイト単位のハミング距離を返す。長さの不一致はinvalid_argument。時間 O(N)、追加領域 O(1)。
inline Int hammingDistance(std::string_view s, std::string_view t) {
    if (s.size() != t.size())
        throw std::invalid_argument("文字列の長さが一致しません");
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
    auto a = reinterpret_cast<const unsigned char *>(s.data());
    auto b = reinterpret_cast<const unsigned char *>(t.data());
#ifndef CPLIB_HAMMING_DISTANCE_SCALAR
    if (__builtin_cpu_supports("popcnt")) {
#ifndef CPLIB_HAMMING_DISTANCE_AVX2
        if (s.size() >= 64 && __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw"))
            return detail::cplib_hamming_avx512(a, b, s.size());
#endif
        if (s.size() >= 32 && __builtin_cpu_supports("avx2"))
            return detail::cplib_hamming_avx2(a, b, s.size());
    }
#endif
#endif
    Int count = 0;
    for (std::size_t i = 0; i < s.size(); ++i)
        count += s[i] != t[i];
    return count;
}
// 文字配列のハミング距離を返す。時間 O(N)、追加領域 O(1)。
inline Int hammingDistance(const std::vector<char> &s, const std::vector<char> &t) {
    return hammingDistance(std::string_view(s.data(), s.size()), std::string_view(t.data(), t.size()));
}
}
