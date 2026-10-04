#pragma once
#include <cplib/common.hpp>
#include <immintrin.h>
#include <stddef.h>

namespace cplib::detail::fenwick_avx2 {

inline ::cplib::Int cplib_fw16_sum(::cplib::Int a, ::cplib::Int b) {
    /* SIMDと同じ桁あふれ動作を、符号付き整数の未定義動作なしで行います。 */
    ::cplib::Int value;
    __builtin_add_overflow(a, b, &value);
    return value;
}

inline ::cplib::Int cplib_fw16_diff(::cplib::Int a, ::cplib::Int b) {
    /* 減算も下位64bitを保持します。 */
    ::cplib::Int value;
    __builtin_sub_overflow(a, b, &value);
    return value;
}

inline void cplib_fw16_build(::cplib::Int *data, const ::cplib::Int *offsets, ::cplib::Int height,
                             ::cplib::Int n) {
    /* 下段を累積和に変換しながら、完全な16要素ブロックの和を上段へ渡します。 */
    for (::cplib::Int h = 0; h < height; ++h, n >>= 4) {
        ::cplib::Int *base = data + offsets[h];
        ::cplib::Int sum = 0;
        for (::cplib::Int i = 0; i <= n; ++i) {
            if (!(i & 15)) {
                if (i)
                    data[offsets[h + 1] + (i >> 4) - 1] = sum;
                sum = 0;
            }
            ::cplib::Int value = base[i];
            base[i] = sum;
            sum = cplib_fw16_sum(sum, value);
        }
    }
}

__attribute__((target("avx2"))) inline void cplib_fw16_add(::cplib::Int *data,
                                                           const ::cplib::Int *offsets,
                                                           ::cplib::Int height, ::cplib::Int index,
                                                           ::cplib::Int delta) {
    /* 各段の16個の累積和を4要素ずつ更新します。 */
    const __m256i value = _mm256_set1_epi64x((long long)delta);
    for (::cplib::Int h = 0; h < height; ++h, index >>= 4) {
        ::cplib::Int *base = data + offsets[h] + (index & ~(::cplib::Int)15);
        const __m256i pos = _mm256_set1_epi64x(index & 15);
        for (int k = 0; k < 16; k += 4) {
            const __m256i mask =
                _mm256_cmpgt_epi64(_mm256_setr_epi64x(k, k + 1, k + 2, k + 3), pos);
            __m256i *p = (__m256i *)(base + k);
            _mm256_storeu_si256(
                p, _mm256_add_epi64(_mm256_loadu_si256(p), _mm256_and_si256(mask, value)));
        }
    }
}

inline ::cplib::Int cplib_fw16_prefix(const ::cplib::Int *data, const ::cplib::Int *offsets,
                                      ::cplib::Int r) {
    /* 各段から1個ずつ累積和を読み出します。 */
    ::cplib::Int sum = 0;
    for (::cplib::Int h = 0; r; ++h, r >>= 4)
        sum = cplib_fw16_sum(sum, data[offsets[h] + r]);
    return sum;
}

inline ::cplib::Int cplib_fw16_get(const ::cplib::Int *data, const ::cplib::Int *offsets,
                                   ::cplib::Int l, ::cplib::Int r) {
    /* 左右の読み出しを並列化し、同じ祖先に達したら終了します。 */
    ::cplib::Int sum = 0;
    for (::cplib::Int h = 0; l != r; ++h, l >>= 4, r >>= 4)
        sum = cplib_fw16_sum(sum, cplib_fw16_diff(data[offsets[h] + r], data[offsets[h] + l]));
    return sum;
}

}
