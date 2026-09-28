#pragma once
#include <immintrin.h>
#include <cstdint>
#include <cstddef>
namespace cplib::detail::bitset_fuse_arithmetic {
// 元の融合カーネルと同じ桁上がり伝播・境界付きロード。
static inline uint64_t cplib_fuse_add64(uint64_t a, uint64_t b, unsigned *carry) {
uint64_t sum = a+b, total = sum+*carry;
*carry = (sum<a) | (total<sum);
return total;
}
__attribute__((target("avx512f"))) static inline __m512i cplib_fuse_add512(__m512i a, __m512i b, unsigned *carry) {
__m512i sum = _mm512_add_epi64(a,b);
unsigned g = _mm512_cmp_epu64_mask(sum,a,_MM_CMPINT_LT);
unsigned p = _mm512_cmpeq_epi64_mask(sum,_mm512_set1_epi64(-1));
unsigned c = (p+(g<<1)+*carry)^p;
*carry = c>>8;
return _mm512_mask_add_epi64(sum,(__mmask8)c,sum,_mm512_set1_epi64(1));
}
static inline uint64_t cplib_fuse_word(const uint64_t *x, size_t n, ptrdiff_t i) {
return i>=0 && (size_t)i<n ? x[i] : 0;
}
static inline uint64_t cplib_fuse_shift0(const uint64_t *x, size_t n, size_t i, size_t k, int left) {
ptrdiff_t j = (ptrdiff_t)i + (left ? -(ptrdiff_t)(k>>6) : (ptrdiff_t)(k>>6));
unsigned b=k&63;
uint64_t a=cplib_fuse_word(x,n,j);
if (!b) return a;
return left ? (a<<b)|(cplib_fuse_word(x,n,j-1)>>(64-b)) : (a>>b)|(cplib_fuse_word(x,n,j+1)<<(64-b));
}
__attribute__((target("avx2"))) static inline __m256i cplib_fuse_load256(const uint64_t* x,size_t n,ptrdiff_t j) {
if(j>=0 && (size_t)j<=n && n-(size_t)j>=4) return _mm256_loadu_si256((const __m256i*)(x+j));
uint64_t a[4]; for(int q=0;q<4;++q) a[q]=cplib_fuse_word(x,n,j+q);
return _mm256_loadu_si256((const __m256i*)a);
}
__attribute__((target("avx2"))) static inline __m256i cplib_fuse_shift256(const uint64_t* x,size_t n,size_t i,size_t k,int left) {
ptrdiff_t j=(ptrdiff_t)i+(left?-(ptrdiff_t)(k>>6):(ptrdiff_t)(k>>6));unsigned b=k&63;
__m256i a=cplib_fuse_load256(x,n,j);if(!b)return a;
__m128i s=_mm_cvtsi32_si128(b),t=_mm_cvtsi32_si128(64-b);
if(left)return _mm256_or_si256(_mm256_sll_epi64(a,s),_mm256_srl_epi64(cplib_fuse_load256(x,n,j-1),t));
return _mm256_or_si256(_mm256_srl_epi64(a,s),_mm256_sll_epi64(cplib_fuse_load256(x,n,j+1),t));
}
__attribute__((target("avx512f"))) static inline __m512i cplib_fuse_load512(const uint64_t* x,size_t n,ptrdiff_t j) {
if(j>=0 && (size_t)j<=n && n-(size_t)j>=8) return _mm512_loadu_si512((const void*)(x+j));
uint64_t a[8]; for(int q=0;q<8;++q) a[q]=cplib_fuse_word(x,n,j+q);
return _mm512_loadu_si512((const void*)a);
}
__attribute__((target("avx512f"))) static inline __m512i cplib_fuse_shift512(const uint64_t* x,size_t n,size_t i,size_t k,int left) {
ptrdiff_t j=(ptrdiff_t)i+(left?-(ptrdiff_t)(k>>6):(ptrdiff_t)(k>>6));unsigned b=k&63;
__m512i a=cplib_fuse_load512(x,n,j);if(!b)return a;
__m128i s=_mm_cvtsi32_si128(b),t=_mm_cvtsi32_si128(64-b);
if(left)return _mm512_or_si512(_mm512_sll_epi64(a,s),_mm512_srl_epi64(cplib_fuse_load512(x,n,j-1),t));
return _mm512_or_si512(_mm512_srl_epi64(a,s),_mm512_sll_epi64(cplib_fuse_load512(x,n,j+1),t));
}
}
