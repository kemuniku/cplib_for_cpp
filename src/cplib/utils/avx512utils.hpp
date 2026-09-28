#pragma once
#include <cplib/common.hpp>
#include <stdexcept>
#include <type_traits>
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__)) && !defined(CPLIB_AVX512UTILS_SCALAR)
#define CPLIB_AVX512UTILS_NATIVE 1

#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
namespace cplib::detail::avx512utils {

inline int cplib_avx512_available(int bytes) {
    // OSのレジスタ保存対応も含め、必要な命令セットを確認する。O(1)。
    return __builtin_cpu_supports("avx512f") &&
        (bytes >= 4 || __builtin_cpu_supports("avx512bw"));
}

#define CPLIB_AVX_BINARY_CASE(T, BITS, OP, EXPR) \
    case OP: \
        for (; i + lanes <= n; i += lanes) { \
            __m512i x = _mm512_loadu_si512((const void*)(a + i * sizeof(T))); \
            __m512i y = _mm512_loadu_si512((const void*)(b + i * sizeof(T))); \
            _mm512_storeu_si512((void*)(d + i * sizeof(T)), EXPR); \
        } \
        if (i < n) { \
            uint64_t mask = UINT64_MAX >> (64 - (n - i)); \
            __m512i x = _mm512_maskz_loadu_epi##BITS(mask, (const void*)(a + i * sizeof(T))); \
            __m512i y = _mm512_maskz_loadu_epi##BITS(mask, (const void*)(b + i * sizeof(T))); \
            _mm512_mask_storeu_epi##BITS((void*)(d + i * sizeof(T)), mask, EXPR); \
        } \
        break;

#define CPLIB_AVX_KERNEL(NAME, T, BITS, SUFFIX, TARGET) \
__attribute__((target(TARGET))) \
inline void cplib_avx512_##NAME(void *aa, void *bb, void *dd, size_t n, int op) { \
    /* 同じ位置の整数演算を行い、端数はマスク付き命令で処理する。O(n)。 */ \
    const unsigned char *a = (const unsigned char*)aa, *b = (const unsigned char*)bb; \
    unsigned char *d = (unsigned char*)dd; \
    const size_t lanes = 64 / sizeof(T); \
    size_t i = 0; \
    switch (op) { \
        CPLIB_AVX_BINARY_CASE(T, BITS, 0, _mm512_add_epi##BITS(x, y)) \
        CPLIB_AVX_BINARY_CASE(T, BITS, 1, _mm512_sub_epi##BITS(x, y)) \
        CPLIB_AVX_BINARY_CASE(T, BITS, 2, _mm512_min_##SUFFIX(x, y)) \
        CPLIB_AVX_BINARY_CASE(T, BITS, 3, _mm512_max_##SUFFIX(x, y)) \
        CPLIB_AVX_BINARY_CASE(T, BITS, 4, _mm512_and_si512(x, y)) \
        CPLIB_AVX_BINARY_CASE(T, BITS, 5, _mm512_or_si512(x, y)) \
        CPLIB_AVX_BINARY_CASE(T, BITS, 6, _mm512_xor_si512(x, y)) \
    } \
} \
__attribute__((target(TARGET))) \
inline void cplib_avx512_reduce_##NAME(void *aa, void *dd, size_t n, int maximum) { \
    /* 複数レーンのmin/maxを集約する。n>0、O(n)。 */ \
    const unsigned char *a = (const unsigned char*)aa; \
    const size_t lanes = 64 / sizeof(T); \
    size_t i = 0; T best; memcpy(&best, a, sizeof(T)); \
    if (n >= lanes) { \
        __m512i v = _mm512_loadu_si512((const void*)a); \
        i = lanes; \
        if (maximum) { \
            for (; i + lanes <= n; i += lanes) \
                v = _mm512_max_##SUFFIX(v, _mm512_loadu_si512((const void*)(a + i * sizeof(T)))); \
        } else { \
            for (; i + lanes <= n; i += lanes) \
                v = _mm512_min_##SUFFIX(v, _mm512_loadu_si512((const void*)(a + i * sizeof(T)))); \
        } \
        T values[64 / sizeof(T)]; _mm512_storeu_si512((void*)values, v); \
        for (size_t j = 0; j < lanes; ++j) \
            if (maximum ? values[j] > best : values[j] < best) best = values[j]; \
    } \
    for (; i < n; ++i) { \
        T x; memcpy(&x, a + i * sizeof(T), sizeof(T)); \
        if (maximum ? x > best : x < best) best = x; \
    } \
    memcpy(dd, &best, sizeof(T)); \
}
CPLIB_AVX_KERNEL(i8, int8_t, 8, epi8, "avx512f,avx512bw")
CPLIB_AVX_KERNEL(u8, uint8_t, 8, epu8, "avx512f,avx512bw")
CPLIB_AVX_KERNEL(i16, int16_t, 16, epi16, "avx512f,avx512bw")
CPLIB_AVX_KERNEL(u16, uint16_t, 16, epu16, "avx512f,avx512bw")
CPLIB_AVX_KERNEL(i32, int32_t, 32, epi32, "avx512f")
CPLIB_AVX_KERNEL(u32, uint32_t, 32, epu32, "avx512f")
CPLIB_AVX_KERNEL(i64, int64_t, 64, epi64, "avx512f")
CPLIB_AVX_KERNEL(u64, uint64_t, 64, epu64, "avx512f")
#undef CPLIB_AVX_KERNEL
#undef CPLIB_AVX_BINARY_CASE
}

#endif
namespace cplib {
template<class T> concept Avx512Integer=std::is_integral_v<T> && !std::is_same_v<T,bool> && sizeof(T)<=8;
template<Avx512Integer T> bool avx512Available(){
#ifdef CPLIB_AVX512UTILS_NATIVE
return detail::avx512utils::cplib_avx512_available(sizeof(T));
#else
return false;
#endif
}
namespace detail::avx512utils {
template<int Op,Avx512Integer T> void binary(std::span<const T> a,std::span<const T> b,std::span<T> dst){if(a.size()!=b.size()||a.size()!=dst.size())throw std::invalid_argument("array lengths differ");if(a.empty())return;
#ifdef CPLIB_AVX512UTILS_NATIVE
if(avx512Available<T>()){
if constexpr(sizeof(T)==1 && std::is_signed_v<T> ==true)cplib_avx512_i8(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
if constexpr(sizeof(T)==1 && std::is_signed_v<T> ==false)cplib_avx512_u8(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
if constexpr(sizeof(T)==2 && std::is_signed_v<T> ==true)cplib_avx512_i16(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
if constexpr(sizeof(T)==2 && std::is_signed_v<T> ==false)cplib_avx512_u16(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
if constexpr(sizeof(T)==4 && std::is_signed_v<T> ==true)cplib_avx512_i32(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
if constexpr(sizeof(T)==4 && std::is_signed_v<T> ==false)cplib_avx512_u32(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
if constexpr(sizeof(T)==8 && std::is_signed_v<T> ==true)cplib_avx512_i64(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
if constexpr(sizeof(T)==8 && std::is_signed_v<T> ==false)cplib_avx512_u64(const_cast<T*>(a.data()),const_cast<T*>(b.data()),dst.data(),a.size(),Op);
return;}
#endif
using U=std::make_unsigned_t<T>;for(std::size_t i=0;i<a.size();++i){T x=a[i],y=b[i];if constexpr(Op==0)dst[i]=std::bit_cast<T>(U(U(x)+U(y)));else if constexpr(Op==1)dst[i]=std::bit_cast<T>(U(U(x)-U(y)));else if constexpr(Op==2)dst[i]=std::min(x,y);else if constexpr(Op==3)dst[i]=std::max(x,y);else if constexpr(Op==4)dst[i]=x&y;else if constexpr(Op==5)dst[i]=x|y;else dst[i]=x^y;}}
template<bool Maximum,Avx512Integer T> T reduce(std::span<const T> a){if(a.empty())throw std::invalid_argument("empty reduction");T out=a[0];
#ifdef CPLIB_AVX512UTILS_NATIVE
if(avx512Available<T>()){
if constexpr(sizeof(T)==1 && std::is_signed_v<T> == true)cplib_avx512_reduce_i8(const_cast<T*>(a.data()),&out,a.size(),Maximum);
if constexpr(sizeof(T)==1 && std::is_signed_v<T> == false)cplib_avx512_reduce_u8(const_cast<T*>(a.data()),&out,a.size(),Maximum);
if constexpr(sizeof(T)==2 && std::is_signed_v<T> == true)cplib_avx512_reduce_i16(const_cast<T*>(a.data()),&out,a.size(),Maximum);
if constexpr(sizeof(T)==2 && std::is_signed_v<T> == false)cplib_avx512_reduce_u16(const_cast<T*>(a.data()),&out,a.size(),Maximum);
if constexpr(sizeof(T)==4 && std::is_signed_v<T> == true)cplib_avx512_reduce_i32(const_cast<T*>(a.data()),&out,a.size(),Maximum);
if constexpr(sizeof(T)==4 && std::is_signed_v<T> == false)cplib_avx512_reduce_u32(const_cast<T*>(a.data()),&out,a.size(),Maximum);
if constexpr(sizeof(T)==8 && std::is_signed_v<T> == true)cplib_avx512_reduce_i64(const_cast<T*>(a.data()),&out,a.size(),Maximum);
if constexpr(sizeof(T)==8 && std::is_signed_v<T> == false)cplib_avx512_reduce_u64(const_cast<T*>(a.data()),&out,a.size(),Maximum);
return out;}
#endif
for(T x:a)if(Maximum?out<x:x<out)out=x;return out;}
}
template<class A,class B,class D> void avx512Add(const A& a,const B& b,D& dst){using T=typename A::value_type;detail::avx512utils::binary<0,T>(a,b,dst);}
template<class A,class B> auto avx512Add(const A& a,const B& b){using T=typename A::value_type;std::vector<T> out(a.size());detail::avx512utils::binary<0,T>(a,b,out);return out;}
template<Avx512Integer T,std::size_t N> auto avx512Add(const std::array<T,N>& a,const std::array<T,N>& b){std::array<T,N> out;detail::avx512utils::binary<0,T>(a,b,out);return out;}
template<class A,class B,class D> void avx512Sub(const A& a,const B& b,D& dst){using T=typename A::value_type;detail::avx512utils::binary<1,T>(a,b,dst);}
template<class A,class B> auto avx512Sub(const A& a,const B& b){using T=typename A::value_type;std::vector<T> out(a.size());detail::avx512utils::binary<1,T>(a,b,out);return out;}
template<Avx512Integer T,std::size_t N> auto avx512Sub(const std::array<T,N>& a,const std::array<T,N>& b){std::array<T,N> out;detail::avx512utils::binary<1,T>(a,b,out);return out;}
template<class A,class B,class D> void avx512Min(const A& a,const B& b,D& dst){using T=typename A::value_type;detail::avx512utils::binary<2,T>(a,b,dst);}
template<class A,class B> auto avx512Min(const A& a,const B& b){using T=typename A::value_type;std::vector<T> out(a.size());detail::avx512utils::binary<2,T>(a,b,out);return out;}
template<Avx512Integer T,std::size_t N> auto avx512Min(const std::array<T,N>& a,const std::array<T,N>& b){std::array<T,N> out;detail::avx512utils::binary<2,T>(a,b,out);return out;}
template<class A,class B,class D> void avx512Max(const A& a,const B& b,D& dst){using T=typename A::value_type;detail::avx512utils::binary<3,T>(a,b,dst);}
template<class A,class B> auto avx512Max(const A& a,const B& b){using T=typename A::value_type;std::vector<T> out(a.size());detail::avx512utils::binary<3,T>(a,b,out);return out;}
template<Avx512Integer T,std::size_t N> auto avx512Max(const std::array<T,N>& a,const std::array<T,N>& b){std::array<T,N> out;detail::avx512utils::binary<3,T>(a,b,out);return out;}
template<class A,class B,class D> void avx512And(const A& a,const B& b,D& dst){using T=typename A::value_type;detail::avx512utils::binary<4,T>(a,b,dst);}
template<class A,class B> auto avx512And(const A& a,const B& b){using T=typename A::value_type;std::vector<T> out(a.size());detail::avx512utils::binary<4,T>(a,b,out);return out;}
template<Avx512Integer T,std::size_t N> auto avx512And(const std::array<T,N>& a,const std::array<T,N>& b){std::array<T,N> out;detail::avx512utils::binary<4,T>(a,b,out);return out;}
template<class A,class B,class D> void avx512Or(const A& a,const B& b,D& dst){using T=typename A::value_type;detail::avx512utils::binary<5,T>(a,b,dst);}
template<class A,class B> auto avx512Or(const A& a,const B& b){using T=typename A::value_type;std::vector<T> out(a.size());detail::avx512utils::binary<5,T>(a,b,out);return out;}
template<Avx512Integer T,std::size_t N> auto avx512Or(const std::array<T,N>& a,const std::array<T,N>& b){std::array<T,N> out;detail::avx512utils::binary<5,T>(a,b,out);return out;}
template<class A,class B,class D> void avx512Xor(const A& a,const B& b,D& dst){using T=typename A::value_type;detail::avx512utils::binary<6,T>(a,b,dst);}
template<class A,class B> auto avx512Xor(const A& a,const B& b){using T=typename A::value_type;std::vector<T> out(a.size());detail::avx512utils::binary<6,T>(a,b,out);return out;}
template<Avx512Integer T,std::size_t N> auto avx512Xor(const std::array<T,N>& a,const std::array<T,N>& b){std::array<T,N> out;detail::avx512utils::binary<6,T>(a,b,out);return out;}
template<class A> auto avx512Min(const A& a){return detail::avx512utils::reduce<false,typename A::value_type>(a);}
template<class A> auto avx512Max(const A& a){return detail::avx512utils::reduce<true,typename A::value_type>(a);}
}
