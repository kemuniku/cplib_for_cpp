#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <cplib/common.hpp>
#include <array>
#include <cstring>
#include <stdexcept>
namespace cplib::detail::combined_collections_native {

#define CPLIB_BS_AVX2 __attribute__((target("avx2")))

#define CPLIB_BS_BINARY(name, scalar, vector) \
CPLIB_BS_AVX2 static inline void name(uint64_t *dst, const uint64_t *x, \
                           const uint64_t *y, size_t n) { \
/* 256ビットずつ論理演算し、残りを64ビットずつ処理します。 */ \
size_t i = 0; \
for (; i + 4 <= n; i += 4) { \
    __m256i a = _mm256_loadu_si256((const __m256i *)(x + i)); \
    [[maybe_unused]] __m256i b = _mm256_loadu_si256((const __m256i *)(y + i)); \
    _mm256_storeu_si256((__m256i *)(dst + i), vector(a, b)); \
} \
for (; i < n; ++i) dst[i] = x[i] scalar y[i]; \
}
CPLIB_BS_BINARY(cplib_bs_and, &, _mm256_and_si256)
CPLIB_BS_BINARY(cplib_bs_or, |, _mm256_or_si256)
CPLIB_BS_BINARY(cplib_bs_xor, ^, _mm256_xor_si256)
#undef CPLIB_BS_BINARY

CPLIB_BS_AVX2 static inline void cplib_bs_not(uint64_t *dst, const uint64_t *x, size_t n) {
/* 256ビットずつ反転します。 */
const __m256i ones = _mm256_set1_epi64x(-1);
size_t i = 0;
for (; i + 4 <= n; i += 4)
    _mm256_storeu_si256((__m256i *)(dst + i), _mm256_xor_si256(
        _mm256_loadu_si256((const __m256i *)(x + i)), ones));
for (; i < n; ++i) dst[i] = ~x[i];
}

CPLIB_BS_AVX2 static inline void cplib_bs_shl(uint64_t *dst, const uint64_t *x,
                                 size_t n, size_t shift) {
/* ゼロ初期化済みの別領域へ左シフトし、隣接ワードからの桁上がりも処理します。 */
const size_t offset = shift >> 6;
const unsigned bits = shift & 63;
const size_t count = n - offset;
size_t i = 0;
if (bits == 0) {
    for (; i + 4 <= count; i += 4)
        _mm256_storeu_si256((__m256i *)(dst + offset + i),
            _mm256_loadu_si256((const __m256i *)(x + i)));
    for (; i < count; ++i) dst[offset + i] = x[i];
    return;
}
const __m128i left = _mm_cvtsi32_si128(bits);
const __m128i right = _mm_cvtsi32_si128(64 - bits);
dst[offset] = x[0] << bits;
i = 1;
for (; i + 4 <= count; i += 4) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(x + i));
    [[maybe_unused]] __m256i b = _mm256_loadu_si256((const __m256i *)(x + i - 1));
    _mm256_storeu_si256((__m256i *)(dst + offset + i), _mm256_or_si256(
        _mm256_sll_epi64(a, left), _mm256_srl_epi64(b, right)));
}
for (; i < count; ++i)
    dst[offset + i] = (x[i] << bits) | (x[i - 1] >> (64 - bits));
}

CPLIB_BS_AVX2 static inline void cplib_bs_shr(uint64_t *dst, const uint64_t *x,
                                 size_t n, size_t shift) {
/* ゼロ初期化済みの別領域へ右シフトし、隣接ワードからの桁下がりも処理します。 */
const size_t offset = shift >> 6;
const unsigned bits = shift & 63;
const size_t count = n - offset;
size_t i = 0;
if (bits == 0) {
    for (; i + 4 <= count; i += 4)
        _mm256_storeu_si256((__m256i *)(dst + i),
            _mm256_loadu_si256((const __m256i *)(x + offset + i)));
    for (; i < count; ++i) dst[i] = x[offset + i];
    return;
}
const __m128i right = _mm_cvtsi32_si128(bits);
const __m128i left = _mm_cvtsi32_si128(64 - bits);
for (; i + 4 < count; i += 4) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(x + offset + i));
    [[maybe_unused]] __m256i b = _mm256_loadu_si256((const __m256i *)(x + offset + i + 1));
    _mm256_storeu_si256((__m256i *)(dst + i), _mm256_or_si256(
        _mm256_srl_epi64(a, right), _mm256_sll_epi64(b, left)));
}
for (; i + 1 < count; ++i)
    dst[i] = (x[offset + i] >> bits) | (x[offset + i + 1] << (64 - bits));
dst[count - 1] = x[n - 1] >> bits;
}

CPLIB_BS_AVX2 static inline __m256i cplib_bs_byte_counts(__m256i x) {
/* 4ビットの参照表から各バイトの立っているビット数を求めます。 */
const __m256i table = _mm256_setr_epi8(
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
const __m256i mask = _mm256_set1_epi8(15);
return _mm256_add_epi8(
    _mm256_shuffle_epi8(table, _mm256_and_si256(x, mask)),
    _mm256_shuffle_epi8(table, _mm256_and_si256(_mm256_srli_epi16(x, 4), mask)));
}

#define CPLIB_BS_COUNT(name, scalar, vector) \
CPLIB_BS_AVX2 static inline size_t name(const uint64_t *x, const uint64_t *y, size_t n) { \
/* 16ベクトルごとにバイトの和を64ビットへ集約し、桁あふれを防ぎます。 */ \
__m256i total = _mm256_setzero_si256(); \
size_t i = 0; \
while (i + 4 <= n) { \
    __m256i local = _mm256_setzero_si256(); \
    size_t end = n - i < 64 ? n : i + 64; \
    for (; i + 4 <= end; i += 4) { \
        __m256i a = _mm256_loadu_si256((const __m256i *)(x + i)); \
        [[maybe_unused]] __m256i b = _mm256_loadu_si256((const __m256i *)(y + i)); \
        local = _mm256_add_epi8(local, cplib_bs_byte_counts(vector)); \
    } \
    total = _mm256_add_epi64(total, _mm256_sad_epu8(local, _mm256_setzero_si256())); \
} \
uint64_t lanes[4]; \
_mm256_storeu_si256((__m256i *)lanes, total); \
size_t result = lanes[0] + lanes[1] + lanes[2] + lanes[3]; \
for (; i < n; ++i) result += __builtin_popcountll(scalar); \
return result; \
}
CPLIB_BS_COUNT(cplib_bs_popcount, x[i], a)
CPLIB_BS_COUNT(cplib_bs_andpopcount, x[i] & y[i], _mm256_and_si256(a, b))
CPLIB_BS_COUNT(cplib_bs_orpopcount, x[i] | y[i], _mm256_or_si256(a, b))
CPLIB_BS_COUNT(cplib_bs_xorpopcount, x[i] ^ y[i], _mm256_xor_si256(a, b))
#undef CPLIB_BS_COUNT

CPLIB_BS_AVX2 static inline uint32_t cplib_bs_bool_mask(const unsigned char *src) {
    /* 32個のboolを比較し、非ゼロの位置を32ビットのマスクに詰めます。 */
    __m256i values = _mm256_loadu_si256((const __m256i *)src);
    return ~(uint32_t)_mm256_movemask_epi8(
        _mm256_cmpeq_epi8(values, _mm256_setzero_si256()));
}

CPLIB_BS_AVX2 static inline void cplib_bs_from_bools(uint64_t *dst, const void *source,
                                            size_t length, size_t words) {
    /* 入力を64ビットずつ詰め、端数と入力より後ろのワードもすべて書き込みます。 */
    const unsigned char *src = (const unsigned char *)source;
    size_t i = 0, word = 0;
    for (; i + 64 <= length; i += 64) {
        uint64_t low = cplib_bs_bool_mask(src + i);
        uint64_t high = cplib_bs_bool_mask(src + i + 32);
        dst[word++] = low | (high << 32);
    }
    if (i < length) {
        uint64_t value = 0;
        size_t j = 0;
        if (length - i >= 32) {
            value = cplib_bs_bool_mask(src + i);
            j = 32;
        }
        for (; j < length - i; ++j)
            value |= (uint64_t)(src[i + j] != 0) << j;
        dst[word++] = value;
    }
    for (; word < words; ++word) dst[word] = 0;
}
#undef CPLIB_BS_AVX2
}
namespace cplib {
// 展開済みファイルに保存された旧固定長AVX2版。検索はワード走査、シフトは別領域へ書き出す。
template<Int N> class BitSet {
 static_assert(N>=0);static constexpr Int Words=(N>>6)+bool(N&63);std::array<UInt,Words> bits_{};
 void check(Int i)const{if(i<0||i>=N)throw std::out_of_range("BitSet index out of bounds");}void trim(){if constexpr(N%64)bits_.back()&=(UInt(1)<<(N%64))-1;}
 inline static constexpr auto ByteStrings=[](){std::array<std::array<char,8>,256> out{};for(Int i=0;i<256;++i)for(Int j=0;j<8;++j)out[i][j]=char('0'+((i>>(7-j))&1));return out;}();
public:
 BitSet()=default;explicit BitSet(std::span<const bool> values){if(Int(values.size())>N)throw std::invalid_argument("initial value is longer than BitSet size");if constexpr(N>0)detail::combined_collections_native::cplib_bs_from_bools(bits_.data(),values.data(),values.size(),Words);}
 explicit BitSet(const std::vector<bool>& values){if(Int(values.size())>N)throw std::invalid_argument("initial value is longer than BitSet size");for(Int i=0;i<Int(values.size());++i)if(values[i])bits_[i>>6]|=UInt(1)<<(i&63);}
 static BitSet fromIndexes(std::span<const Int> indexes){BitSet out;for(Int i:indexes){out.check(i);out.bits_[i>>6]|=UInt(1)<<(i&63);}return out;}
 static constexpr Int len(){return N;}bool operator[](Int i)const{check(i);return (bits_[i>>6]>>(i&63))&1;}void set(Int i,bool value){check(i);if(value)bits_[i>>6]|=UInt(1)<<(i&63);else bits_[i>>6]&=~(UInt(1)<<(i&63));}
 struct Reference {BitSet* owner;Int i;operator bool()const{return std::as_const(*owner)[i];}Reference& operator=(bool value){owner->set(i,value);return *this;}Reference& operator=(Int value){if(value==0||value==1)owner->set(i,value==1);return *this;}Reference& operator=(int value){return *this=Int(value);}Reference& operator=(const Reference& value){return *this=bool(value);}};Reference operator[](Int i){return {this,i};}
 BitSet& operator&=(const BitSet& b){if constexpr(N>0)detail::combined_collections_native::cplib_bs_and(bits_.data(),bits_.data(),b.bits_.data(),Words);return *this;}BitSet& operator|=(const BitSet& b){if constexpr(N>0)detail::combined_collections_native::cplib_bs_or(bits_.data(),bits_.data(),b.bits_.data(),Words);return *this;}BitSet& operator^=(const BitSet& b){if constexpr(N>0)detail::combined_collections_native::cplib_bs_xor(bits_.data(),bits_.data(),b.bits_.data(),Words);return *this;}
 friend BitSet operator&(const BitSet& a,const BitSet& b){BitSet out;if constexpr(N>0)detail::combined_collections_native::cplib_bs_and(out.bits_.data(),a.bits_.data(),b.bits_.data(),Words);return out;}friend BitSet operator|(const BitSet& a,const BitSet& b){BitSet out;if constexpr(N>0)detail::combined_collections_native::cplib_bs_or(out.bits_.data(),a.bits_.data(),b.bits_.data(),Words);return out;}friend BitSet operator^(const BitSet& a,const BitSet& b){BitSet out;if constexpr(N>0)detail::combined_collections_native::cplib_bs_xor(out.bits_.data(),a.bits_.data(),b.bits_.data(),Words);return out;}
 BitSet operator<<(Int count)const{if(count<0)throw std::invalid_argument("shift count must be non-negative");BitSet out;if constexpr(N>0)if(count<N){detail::combined_collections_native::cplib_bs_shl(out.bits_.data(),bits_.data(),Words,count);out.trim();}return out;}BitSet operator>>(Int count)const{if(count<0)throw std::invalid_argument("shift count must be non-negative");BitSet out;if constexpr(N>0)if(count<N)detail::combined_collections_native::cplib_bs_shr(out.bits_.data(),bits_.data(),Words,count);return out;}
 BitSet operator~()const{BitSet out;if constexpr(N>0){detail::combined_collections_native::cplib_bs_not(out.bits_.data(),bits_.data(),Words);out.trim();}return out;}
 Int popcount()const{if constexpr(N>0)return detail::combined_collections_native::cplib_bs_popcount(bits_.data(),bits_.data(),Words);else return 0;}Int andpopcount(const BitSet& b)const{if constexpr(N>0)return detail::combined_collections_native::cplib_bs_andpopcount(bits_.data(),b.bits_.data(),Words);else return 0;}Int orpopcount(const BitSet& b)const{if constexpr(N>0)return detail::combined_collections_native::cplib_bs_orpopcount(bits_.data(),b.bits_.data(),Words);else return 0;}Int xorpopcount(const BitSet& b)const{if constexpr(N>0)return detail::combined_collections_native::cplib_bs_xorpopcount(bits_.data(),b.bits_.data(),Words);else return 0;}
 struct Iterator {const BitSet* owner;Int word;UInt remaining;void skip(){while(!remaining&&word<Words){++word;if(word<Words)remaining=owner->bits_[word];}}Int operator*()const{return word*64+std::countr_zero(remaining);}Iterator& operator++(){remaining&=remaining-1;skip();return *this;}bool operator!=(const Iterator& b)const{return word!=b.word||remaining!=b.remaining;}};
 Iterator begin()const{Iterator out{this,0,Words?bits_[0]:0};out.skip();return out;}Iterator end()const{return {this,Words,0};}const BitSet& items()const{return *this;}
 Int lowestBit()const{for(Int i=0;i<Words;++i)if(bits_[i])return i*64+std::countr_zero(bits_[i]);return -1;}
 std::string str()const{std::string out(N,'\0');Int finish=N;for(Int i=0;i<Words;++i){UInt word=bits_[i];unsigned remaining=unsigned(std::min(Int(64),N-i*64));while(remaining>=8){finish-=8;std::memcpy(out.data()+finish,ByteStrings[word&255].data(),8);word>>=8;remaining-=8;}for(unsigned j=0;j<remaining;++j){out[--finish]=char('0'+(word&1));word>>=1;}}return out;}
};
template<Int N> auto initBitSet(){return BitSet<N>();}template<Int N,class R> auto initBitSet(const R& values){return BitSet<N>(values);}template<Int N> auto initBitSetFromIndexes(std::span<const Int> indexes){return BitSet<N>::fromIndexes(indexes);}template<Int N> Int len(const BitSet<N>& a){return a.len();}template<Int N> Int popcount(const BitSet<N>& a){return a.popcount();}template<Int N> Int andpopcount(const BitSet<N>& a,const BitSet<N>& b){return a.andpopcount(b);}template<Int N> Int orpopcount(const BitSet<N>& a,const BitSet<N>& b){return a.orpopcount(b);}template<Int N> Int xorpopcount(const BitSet<N>& a,const BitSet<N>& b){return a.xorpopcount(b);}template<Int N> Int lowestBit(const BitSet<N>& a){return a.lowestBit();}
}
