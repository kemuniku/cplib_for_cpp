#pragma once
#include <cplib/math/inner_math.hpp>
#include <cplib/math/isprime.hpp>
#include <cplib/str/run_length_encode.hpp>
#include <random>
#include <unordered_map>
namespace cplib {
namespace detail {
inline thread_local std::mt19937_64 factor_rng{std::random_device{}()};
// Brent版Pollard rhoで素因数を探す。期待O(p^(1/2))、pは最小素因数。
inline Int find_factor(Int n) {
    if (!(n&1)) return 2;
    if (isprime(n)) return n;
    constexpr Int m=128;
    for (;;) {
        Int x=1,ys=1,q=1,r=1,g=1;
        Int rnd=std::uniform_int_distribution<Int>(2,n-1)(factor_rng);
        Int y=std::uniform_int_distribution<Int>(2,n-1)(factor_rng);
        auto f=[&](Int v){return add(mul(v,v,n),rnd,n);};
        while (g==1) {
            x=y;
            for (Int i=0;i<r;++i) y=f(y);
            for (Int k=0;k<r;k+=m) {
                ys=y;
                for (Int i=0;i<std::min(m,r-k);++i) { y=f(y);q=mul(q,std::abs(x-y),n); }
                g=std::gcd(q,n); if(g!=1) break;
            }
            r*=2;
        }
        if (g==n) {g=1;while(g==1){ys=f(ys);g=std::gcd(n,std::abs(x-ys));}}
        if (g<n) {
            if (isprime(g)) return g;
            if (isprime(n/g)) return n/g;
            return find_factor(g);
        }
    }
}
}
// 素因数を重複込みで列挙する。元実装と同じPollard rho法。
inline std::vector<Int> primefactor(Int n, bool sorted=true) {
    std::vector<Int> result;
    while(n>1 && !isprime(n)) {
        Int p=detail::find_factor(n);
        while(n%p==0){result.push_back(p);n/=p;}
    }
    if(n>1) result.push_back(n);
    if(sorted) std::sort(result.begin(),result.end());
    return result;
}
// 素因数ごとの指数を返す。
inline std::unordered_map<Int,Int> primefactor_table(Int n) {
    std::unordered_map<Int,Int> result;
    for(auto p:primefactor(n)) ++result[p];
    return result;
}
// 素因数と指数の組を昇順で返す。
inline std::vector<std::pair<Int,Int>> primefactor_tuple(Int n) {return run_length_encode(primefactor(n));}
}
