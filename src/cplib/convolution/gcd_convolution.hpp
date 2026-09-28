#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 素数ごとの約数変換。gcd(0,i)=i。O(N log log N)、領域O(N)。
template<class T> std::vector<T> gcdConvolution(const std::vector<T>& a,const std::vector<T>& b){assert(a.size()==b.size());Int n=a.size();if(!n)return {};auto out=a,right=b;std::vector<bool> composite(n);std::vector<Int> primes;for(Int p=2;p<n;++p)if(!composite[p]){primes.push_back(p);for(Int i=(n-1)/p;i>=1;--i){composite[i*p]=true;out[i]+=out[i*p];right[i]+=right[i*p];}}for(Int i=0;i<n;++i)out[i]*=right[i];for(Int p:primes)for(Int i=1;i<=(n-1)/p;++i)out[i]-=out[i*p];for(Int i=1;i<n;++i)out[i]+=a[0]*b[i]+a[i]*b[0];return out;}
}
