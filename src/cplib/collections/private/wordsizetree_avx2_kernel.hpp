#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
namespace cplib::detail::wordsize_avx2 {




#define WST_AVX __attribute__((target("avx2")))
/* 各ノードの256ビットを、連続する4個の64ビット整数に格納します。 */
WST_AVX static inline unsigned wst_lanes(const uint64_t *p) {
    /* 非零の64ビット整数の位置をビットマスクで返します。 */
    __m256i v = _mm256_loadu_si256((const __m256i *)p);
    return (~(unsigned)_mm256_movemask_pd(_mm256_castsi256_pd(
        _mm256_cmpeq_epi64(v, _mm256_setzero_si256())))) & 15u;
}
WST_AVX static inline int wst_next(const uint64_t *p, int bit) {
    /* ノード内でbit以上の最小の要素を返し、存在しなければ-1を返します。 */
    if (bit >= 256) return -1;
    unsigned lane = (unsigned)bit >> 6;
    uint64_t word = p[lane] & (UINT64_MAX << (bit & 63));
    if (word) return (int)(lane * 64 + __builtin_ctzll(word));
    unsigned mask = wst_lanes(p) & (15u << (lane + 1));
    if (!mask) return -1;
    lane = __builtin_ctz(mask);
    return (int)(lane * 64 + __builtin_ctzll(p[lane]));
}
WST_AVX static inline int wst_prev(const uint64_t *p, int bit) {
    /* ノード内でbit以下の最大の要素を返し、存在しなければ-1を返します。 */
    if (bit < 0) return -1;
    unsigned lane = (unsigned)bit >> 6;
    uint64_t word = p[lane] & (UINT64_MAX >> (63 - (bit & 63)));
    if (word) return (int)(lane * 64 + 63 - __builtin_clzll(word));
    unsigned mask = wst_lanes(p) & ((1u << lane) - 1);
    if (!mask) return -1;
    lane = 31 - __builtin_clz(mask);
    return (int)(lane * 64 + 63 - __builtin_clzll(p[lane]));
}
WST_AVX static inline void wst_init(const void *input, size_t n, unsigned char one,
                            uint64_t *leaf, uint64_t *mid, uint64_t *top) {
    /* oneと一致する位置から、ゼロ初期化された各段のビット集合を構築します。 */
    const unsigned char *v = (const unsigned char *)input;
    size_t i = 0;
    const __m256i target = _mm256_set1_epi8((char)one);
    for (; i + 64 <= n; i += 64) {
        __m256i a = _mm256_loadu_si256((const __m256i *)(v + i));
        __m256i b = _mm256_loadu_si256((const __m256i *)(v + i + 32));
        uint32_t lo = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, target));
        uint32_t hi = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(b, target));
        leaf[i >> 6] = (uint64_t)lo | ((uint64_t)hi << 32);
    }
    if (i < n) {
        uint64_t word = 0;
        for (size_t j = 0; i + j < n; ++j) word |= (uint64_t)(v[i+j] == one) << j;
        leaf[i >> 6] = word;
    }
    size_t nodes = (n + 255) >> 8;
    for (size_t j = 0; j < nodes; ++j)
        mid[j >> 6] |= (uint64_t)(wst_lanes(leaf + j * 4) != 0) << (j & 63);
    for (size_t j = 0; j < ((nodes + 255) >> 8); ++j)
        top[j >> 6] |= (uint64_t)(wst_lanes(mid + j * 4) != 0) << (j & 63);
}
WST_AVX static inline void wst_incl(uint64_t *leaf, uint64_t *mid, uint64_t *top, unsigned x) {
    /* 要素xを追加し、上位のビット集合を更新します。 */
    leaf[x >> 6] |= UINT64_C(1) << (x & 63);
    x >>= 8;
    mid[x >> 6] |= UINT64_C(1) << (x & 63);
    x >>= 8;
    top[x >> 6] |= UINT64_C(1) << (x & 63);
}
WST_AVX static inline void wst_excl(uint64_t *leaf, uint64_t *mid, uint64_t *top, unsigned x) {
    /* 要素xを削除し、空になったノードを上位のビット集合から除きます。 */
    leaf[x >> 6] &= ~(UINT64_C(1) << (x & 63));
    x >>= 8;
    if (wst_lanes(leaf + x * 4)) return;
    mid[x >> 6] &= ~(UINT64_C(1) << (x & 63));
    x >>= 8;
    if (wst_lanes(mid + x * 4)) return;
    top[x >> 6] &= ~(UINT64_C(1) << (x & 63));
}
WST_AVX static inline int wst_ge(const uint64_t *leaf, const uint64_t *mid,
                          const uint64_t *top, unsigned x) {
    /* x以上の最小の要素を返し、存在しなければ-1を返します。 */
    unsigned node = x >> 8;
    int bit = wst_next(leaf + node * 4, x & 255);
    if (bit >= 0) return (int)(node * 256 + bit);
    unsigned parent = node >> 8;
    bit = wst_next(mid + parent * 4, (node & 255) + 1);
    if (bit < 0) {
        int upper = wst_next(top, parent + 1);
        if (upper < 0) return -1;
        parent = (unsigned)upper;
        bit = wst_next(mid + parent * 4, 0);
    }
    node = parent * 256 + bit;
    return (int)(node * 256 + wst_next(leaf + node * 4, 0));
}
WST_AVX static inline int wst_le(const uint64_t *leaf, const uint64_t *mid,
                          const uint64_t *top, unsigned x) {
    /* x以下の最大の要素を返し、存在しなければ-1を返します。 */
    unsigned node = x >> 8;
    int bit = wst_prev(leaf + node * 4, x & 255);
    if (bit >= 0) return (int)(node * 256 + bit);
    unsigned parent = node >> 8;
    bit = wst_prev(mid + parent * 4, (int)(node & 255) - 1);
    if (bit < 0) {
        int upper = wst_prev(top, (int)parent - 1);
        if (upper < 0) return -1;
        parent = (unsigned)upper;
        bit = wst_prev(mid + parent * 4, 255);
    }
    node = parent * 256 + bit;
    return (int)(node * 256 + wst_prev(leaf + node * 4, 255));
}
#undef WST_AVX

}
