#pragma once
#include <cplib/fps/formal_power_series.hpp>
namespace cplib {
// f(x+c)を階乗で重み付けした一回の畳み込みから求める。O(M(N))。
template<Modint T> std::vector<T> taylorShift(const std::vector<T>& f,std::type_identity_t<T> c){Int n=f.size();if(!n)return {};assert(n<=T::umod());std::vector<T> fact(n),factInv(n),left(n),right(n);fact[0]=1;for(Int i=1;i<n;++i)fact[i]=fact[i-1]*T(i);factInv.back()=fact.back().inv();for(Int i=n-1;i>=1;--i)factInv[i-1]=factInv[i]*T(i);T cpow=1;for(Int i=0;i<n;++i){left[n-1-i]=f[i]*fact[i];right[i]=cpow*factInv[i];cpow*=c;}auto product=left*right;std::vector<T> out(n);for(Int i=0;i<n;++i)out[i]=product[n-1-i]*factInv[i];return out;}
}
