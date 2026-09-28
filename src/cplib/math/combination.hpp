#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 階乗・逆元・逆階乗をO(max_N)で前計算し、組合せ等をO(1)で返す。
template<class ModInt> struct Combination_Type {
    std::vector<ModInt> fact,inv,fact_inv;
    explicit Combination_Type(Int max_N):fact(max_N+1),inv(max_N+1),fact_inv(max_N+1){assert(max_N>=0);fact[0]=1;fact_inv[0]=1;if(max_N>=1){fact[1]=1;inv[1]=1;fact_inv[1]=1;}for(Int i=2;i<=max_N;++i){fact[i]=fact[i-1]*i;inv[i]=-inv[Int(ModInt::umod())%i]*(Int(ModInt::umod())/i);fact_inv[i]=fact_inv[i-1]*inv[i];}}
    ModInt ncr(Int n,Int r)const{if(n<0||r<0||n<r)return 0;return fact[n]*fact_inv[n-r]*fact_inv[r];}
    ModInt npr(Int n,Int r)const{if(n<0||r<0||n<r)return 0;return fact[n]*fact_inv[n-r];}
    ModInt nhr(Int n,Int r)const{if(n==0&&r==0)return 1;return ncr(n+r-1,r);}
};
template<class M> auto initCombination(Int max_N){return Combination_Type<M>(max_N);}
template<class M> auto ncr(const Combination_Type<M>& c,Int n,Int r){return c.ncr(n,r);}
template<class M> auto npr(const Combination_Type<M>& c,Int n,Int r){return c.npr(n,r);}
template<class M> auto nhr(const Combination_Type<M>& c,Int n,Int r){return c.nhr(n,r);}
}
