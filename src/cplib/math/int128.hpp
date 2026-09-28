#pragma once
#include <cplib/common.hpp>
#include <cstring>
#include <iostream>
#include <cctype>
#include <cstdio>
#include <functional>
namespace cplib {
namespace detail {


    inline unsigned long long parseuint_raw8b(const unsigned long long &x) {
        // https://zenn.dev/mizar/articles/fc87d667153080
        unsigned long long result = (x & 0x0f0f0f0f0f0f0f0f);
        result *= ((10ul << 8) + 1); result >>= 8; result &= 0x00ff00ff00ff00ff;
        result *= ((100ul << 16) + 1); result >>= 16; result &= 0x0000ffff0000ffff;
        result *= ((10000ul << 32) + 1); result >>= 32;
        return result;
    }
    inline unsigned long long parseint_raw8b_wrap(const char* p, size_t sz) {
        char c[8] = {0};
        for (size_t i=0; i<sz; i++) c[8-sz+i] = *(p++);
        unsigned long long x;
        std::memcpy(&x, &c, 8);
        return parseuint_raw8b(x);
    }
    inline __int128_t parse_int128(const char* p) {
        // 符号付き128ビット整数を負の値として累積し、最小値も安全に解析する。
        bool minus = *p == '-' ? (p++, true) : false;
        const __int128_t base[9] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
        __int128_t result = 0;
        while (1) {
            size_t sz = 0;
            for (size_t i=0; i<8; i++) {
                if (*(p + sz) == '\0') break;
                sz++;
            }
            result = result * base[sz] - parseint_raw8b_wrap(p, sz);
            p += sz;
            if (*p == '\0') break;
        }
        return minus ? result : -result;
    }
    constexpr size_t INT128_DIGIT_STRING_SIZE = 10000;
    constexpr size_t INT128_DIGIT_STRING_LENGHT = 4;
    struct Int128FourDigitStrings {
        char d[INT128_DIGIT_STRING_SIZE * INT128_DIGIT_STRING_LENGHT];
        constexpr Int128FourDigitStrings() : d() {
            for (size_t i=0; i<INT128_DIGIT_STRING_SIZE * INT128_DIGIT_STRING_LENGHT; i++) d[i] = '0';
            for (size_t i=0; i<INT128_DIGIT_STRING_SIZE; i++) {
                size_t pos = INT128_DIGIT_STRING_LENGHT - 1;
                size_t tmp = i;
                while (tmp) {
                    d[i*INT128_DIGIT_STRING_LENGHT+pos--] = "0123456789"[tmp % 10];
                    tmp /= 10;
                }
            }
        }
    };
    inline thread_local char int128_string_buffer[41];
    constexpr auto int128_four_digit_strings = Int128FourDigitStrings();
    inline char* to_string(__int128_t &x) {
        // 絶対値を符号なしで求め、NUL終端付きの十進文字列を返す。
        __uint128_t tmp = static_cast<__uint128_t>(x);
        if (x < 0) tmp = -tmp;
        char* end = std::end(int128_string_buffer) - 1;
        *end = '\0';
        char* d = end;
        while (tmp >= INT128_DIGIT_STRING_SIZE) {
            size_t pos = (tmp % INT128_DIGIT_STRING_SIZE) * INT128_DIGIT_STRING_LENGHT;
            d -= INT128_DIGIT_STRING_LENGHT;
            std::memcpy(d, int128_four_digit_strings.d+pos, INT128_DIGIT_STRING_LENGHT);
            tmp /= INT128_DIGIT_STRING_SIZE;
        }
        while (tmp > 0) {
            *(--d) = "0123456789"[tmp % 10];
            tmp /= 10;
        }
        if (d == end) *(--d) = '0';
        if (x < 0) *(--d) = '-';
        return d;
    }
    inline std::ostream &operator<<(std::ostream &dest, __int128_t &x) {
        // NUL終端付きの文字列をストリームへ出力する。
        return dest << to_string(x);
    }
    inline __int128_t read_and_parse_int128(int) {
        // 空白を読み飛ばして整数を読む。値のないEOFでは0を返す。計算量O(文字数)。
        int c = getchar_unlocked();
        while (c != EOF && std::isspace(static_cast<unsigned char>(c))) {
            c = getchar_unlocked();
        }
        bool minus = c == '-';
        if (c == '-' || c == '+') c = getchar_unlocked();
        __int128_t result = 0;
        while (c >= '0' && c <= '9') {
            result = result * 10 - (c - '0');
            c = getchar_unlocked();
        }
        return minus ? result : -result;
    }
    inline void output_int128(__int128_t &x) { std::cout << x << '\n'; }
    
}
using Int128=__int128_t;
using UInt128=__uint128_t;
// 128ビット整数に変換する。O(1)。
template<class T> constexpr Int128 to_Int128(T x){return static_cast<Int128>(x);}
inline Int to_int(Int128 x){return static_cast<Int>(x);}
// 十進文字列を解析する。O(文字数)。
inline Int128 parseInt128(const std::string& s){return detail::parse_int128(s.c_str());}
inline Int128 read_and_parse_int128(){return detail::read_and_parse_int128(0);}
inline std::string to_string(Int128 x){return detail::to_string(x);}
inline void put(Int128 x){detail::output_int128(x);}
inline Int128 abs(Int128 x){return x>=0?x:-x;}
inline int cmp(Int128 x,Int128 y){return x<y?-1:x==y?0:1;}
// 二分累乗で整数乗を求める。O(log n)。
inline Int128 pow(Int128 x,Int128 n){Int128 result=1;while(n>0){if(n&1)result*=x;if(n>1)x*=x;n>>=1;}return result;}
inline Int128 pow(Int128 x,Int128 n,Int128 m){assert(m!=0);if(m==1)return 0;Int128 result=1;x%=m;while(n>0){if(n&1)result=result*x%m;if(n>1)x=x*x%m;n>>=1;}return result;}
inline std::size_t hash(Int128 x){auto lo=std::hash<UInt>{}(static_cast<UInt>(x));auto hi=std::hash<UInt>{}(static_cast<UInt>(static_cast<UInt128>(x)>>64));return lo^(hi+0x9e3779b97f4a7c15ULL+(lo<<6)+(lo>>2));}
}
