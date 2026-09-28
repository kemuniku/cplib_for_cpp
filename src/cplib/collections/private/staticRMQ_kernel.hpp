#pragma once
#include <cstdint>

#include <immintrin.h>
#include <algorithm>
#include <limits>

namespace cplib::detail::static_rmq_kernel {
constexpr std::int64_t shift = 6;
constexpr std::int64_t block_size = std::int64_t(1) << shift;
#define CPLIB_RMQ_AVX __attribute__((target("avx2")))

static inline bool avx_available() {
    // AVX2対応CPUか判定し、検証時にはスカラー処理を強制できる。
#ifdef CPLIB_STATIC_RMQ_NO_AVX2
    return false;
#else
    return __builtin_cpu_supports("avx2");
#endif
}

template<int Width> struct Lanes;
template<> struct Lanes<4> {
    CPLIB_RMQ_AVX static inline __m256i minimum(__m256i a, __m256i b) {
        // 符号付き32ビット整数を8個まとめて比較する。
        return _mm256_min_epi32(a, b);
    }
    CPLIB_RMQ_AVX static inline __m256i broadcast(std::int32_t x) {
        // 同じ値を全レーンに複製する。
        return _mm256_set1_epi32(x);
    }
    CPLIB_RMQ_AVX static inline __m256i reverse(__m256i x) {
        // レーンの順序を反転する。
        return _mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(7,6,5,4,3,2,1,0));
    }
    CPLIB_RMQ_AVX static inline __m256i prefix(__m256i x) {
        // 8レーン内の累積最小値を並列に計算する。
        const __m256i inf = broadcast(std::numeric_limits<std::int32_t>::max());
        __m256i y = _mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(0,0,1,2,3,4,5,6));
        x = minimum(x, _mm256_blend_epi32(inf, y, 0xfe));
        y = _mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(0,0,0,1,2,3,4,5));
        x = minimum(x, _mm256_blend_epi32(inf, y, 0xfc));
        y = _mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(0,0,0,0,0,1,2,3));
        return minimum(x, _mm256_blend_epi32(inf, y, 0xf0));
    }
    CPLIB_RMQ_AVX static inline __m256i last(__m256i x) {
        // 最後のレーンを全レーンに複製する。
        return _mm256_permutevar8x32_epi32(x, _mm256_set1_epi32(7));
    }
    CPLIB_RMQ_AVX static inline std::int32_t reduce(__m256i x) {
        // 8レーンの最小値を1個の整数にまとめる。
        __m128i y = _mm_min_epi32(_mm256_castsi256_si128(x), _mm256_extracti128_si256(x, 1));
        y = _mm_min_epi32(y, _mm_shuffle_epi32(y, 0x4e));
        y = _mm_min_epi32(y, _mm_shuffle_epi32(y, 0xb1));
        return _mm_cvtsi128_si32(y);
    }
};
template<> struct Lanes<8> {
    CPLIB_RMQ_AVX static inline __m256i minimum(__m256i a, __m256i b) {
        // AVX2にない符号付き64ビット最小値を比較と選択で求める。
        return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(a, b));
    }
    CPLIB_RMQ_AVX static inline __m256i broadcast(std::int64_t x) {
        // 同じ値を全レーンに複製する。
        return _mm256_set1_epi64x(x);
    }
    CPLIB_RMQ_AVX static inline __m256i reverse(__m256i x) {
        // レーンの順序を反転する。
        return _mm256_permute4x64_epi64(x, 0x1b);
    }
    CPLIB_RMQ_AVX static inline __m256i prefix(__m256i x) {
        // 4レーン内の累積最小値を並列に計算する。
        const __m256i inf = broadcast(std::numeric_limits<std::int64_t>::max());
        __m256i y = _mm256_permute4x64_epi64(x, 0x90);
        x = minimum(x, _mm256_blend_epi32(inf, y, 0xfc));
        y = _mm256_permute4x64_epi64(x, 0x40);
        return minimum(x, _mm256_blend_epi32(inf, y, 0xf0));
    }
    CPLIB_RMQ_AVX static inline __m256i last(__m256i x) {
        // 最後のレーンを全レーンに複製する。
        return _mm256_permute4x64_epi64(x, 0xff);
    }
    CPLIB_RMQ_AVX static inline std::int64_t reduce(__m256i x) {
        // 4レーンの最小値を1個の整数にまとめる。
        x = minimum(x, _mm256_permute4x64_epi64(x, 0x4e));
        x = minimum(x, _mm256_permute4x64_epi64(x, 0xb1));
        return _mm_cvtsi128_si64(_mm256_castsi256_si128(x));
    }
};

template<class T> CPLIB_RMQ_AVX static T scan(const T* p, std::int64_t n) {
    // 区間内だけをロードし、端数は末尾のベクトルと重ねて処理する。
    using V = Lanes<sizeof(T)>;
    constexpr std::int64_t width = 32 / sizeof(T);
    __m256i x = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
    for (std::int64_t i = width; i + width <= n; i += width)
        x = V::minimum(x, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + i)));
    if (n % width != 0)
        x = V::minimum(x, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + n - width)));
    return V::reduce(x);
}

template<class T> static void build_scalar(const T* data, T* prefix, T* suffix, T** rows, std::int64_t n) {
    // 各ブロックの累積最小値とブロック間のスパーステーブルを構築する。
    T* table = rows[0];
    const std::int64_t blocks = ((n - 1) >> shift) + 1;
    for (std::int64_t first = 0; first < n; first += block_size) {
        const std::int64_t end = std::min(first + block_size, n);
        prefix[first] = data[first];
        for (std::int64_t i = first + 1; i < end; ++i) prefix[i] = std::min(prefix[i-1], data[i]);
        suffix[end-1] = data[end-1];
        for (std::int64_t i = end - 1; i > first; --i) suffix[i-1] = std::min(suffix[i], data[i-1]);
        table[first >> shift] = prefix[end-1];
    }
    T* previous = table;
    std::int64_t level = 1;
    for (std::int64_t distance = 1; distance * 2 <= blocks; distance *= 2, ++level) {
        T* next = rows[level];
        const std::int64_t length = blocks - distance * 2 + 1;
        for (std::int64_t i = 0; i < length; ++i) next[i] = std::min(previous[i], previous[i + distance]);
        previous = next;
    }
}

template<class T> CPLIB_RMQ_AVX static void build_avx(const T* data, T* prefix, T* suffix, T** rows, std::int64_t n) {
    // 累積最小値とスパーステーブルの構築をAVX2で並列化する。
    using V = Lanes<sizeof(T)>;
    constexpr std::int64_t width = 32 / sizeof(T);
    T* table = rows[0];
    const std::int64_t blocks = ((n - 1) >> shift) + 1;
    std::int64_t first = 0;
    for (; first + block_size <= n; first += block_size) {
        __m256i carry = V::broadcast(std::numeric_limits<T>::max());
        for (std::int64_t i = first; i < first + block_size; i += width) {
            __m256i x = V::prefix(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + i)));
            x = V::minimum(x, carry);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(prefix + i), x);
            carry = V::last(x);
        }
        carry = V::broadcast(std::numeric_limits<T>::max());
        for (std::int64_t i = first + block_size; i > first; i -= width) {
            __m256i x = V::reverse(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + i - width)));
            x = V::minimum(V::prefix(x), carry);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(suffix + i - width), V::reverse(x));
            carry = V::last(x);
        }
        table[first >> shift] = prefix[first + block_size - 1];
    }
    if (first < n) {
        prefix[first] = data[first];
        for (std::int64_t i = first + 1; i < n; ++i) prefix[i] = std::min(prefix[i-1], data[i]);
        suffix[n-1] = data[n-1];
        for (std::int64_t i = n - 1; i > first; --i) suffix[i-1] = std::min(suffix[i], data[i-1]);
        table[first >> shift] = prefix[n-1];
    }
    T* previous = table;
    std::int64_t level = 1;
    for (std::int64_t distance = 1; distance * 2 <= blocks; distance *= 2, ++level) {
        T* next = rows[level];
        const std::int64_t length = blocks - distance * 2 + 1;
        std::int64_t i = 0;
        for (; i + width <= length; i += width) {
            __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(previous + i));
            __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(previous + i + distance));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(next + i), V::minimum(a, b));
        }
        for (; i < length; ++i) next[i] = std::min(previous[i], previous[i + distance]);
        previous = next;
    }
}

template<class T> static void build(const T* data, T* prefix, T* suffix, T** rows, std::int64_t n) {
    // CPUに対応した構築処理を選ぶ。
    if (avx_available()) build_avx(data, prefix, suffix, rows, n);
    else build_scalar(data, prefix, suffix, rows, n);
}

template<class T> static T scan_dispatch(const T* data, std::int64_t n) {
    // CPUに対応したブロック内の最小値走査を選ぶ。
    if (avx_available()) return scan(data, n);
    T result = data[0];
    for (std::int64_t i = 1; i < n; ++i) result = std::min(result, data[i]);
    return result;
}
#undef CPLIB_RMQ_AVX
}

inline void cplib_static_rmq_build_i32(void* data, void* prefix, void* suffix, void* rows, std::int64_t n) {
    // 32ビット整数用の構築処理をNimから呼び出す。
    cplib::detail::static_rmq_kernel::build(static_cast<const std::int32_t*>(data), static_cast<std::int32_t*>(prefix),
        static_cast<std::int32_t*>(suffix), static_cast<std::int32_t**>(rows), n);
}
inline void cplib_static_rmq_build_i64(void* data, void* prefix, void* suffix, void* rows, std::int64_t n) {
    // 64ビット整数用の構築処理をNimから呼び出す。
    cplib::detail::static_rmq_kernel::build(static_cast<const std::int64_t*>(data), static_cast<std::int64_t*>(prefix),
        static_cast<std::int64_t*>(suffix), static_cast<std::int64_t**>(rows), n);
}
inline std::int32_t cplib_static_rmq_scan_i32(void* data, std::int64_t n) {
    // 32ビット整数用の問い合わせをNimから呼び出す。
    return cplib::detail::static_rmq_kernel::scan_dispatch(static_cast<const std::int32_t*>(data), n);
}
inline std::int64_t cplib_static_rmq_scan_i64(void* data, std::int64_t n) {
    // 64ビット整数用の問い合わせをNimから呼び出す。
    return cplib::detail::static_rmq_kernel::scan_dispatch(static_cast<const std::int64_t*>(data), n);
}
