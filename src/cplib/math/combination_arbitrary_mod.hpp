#pragma once
#include <cplib/common.hpp>
namespace cplib {
class CombinationArbitraryMod {
 struct Factor{Int p,q=0,modulus=1,coefficient;std::vector<std::int32_t> fact,inv;};
 Int modulus;std::vector<Factor> factors;
 static Int power(Int a,Int e,Int m){Int r=1;while(e){if(e&1)r=r*a%m;a=a*a%m;e>>=1;}return r;}
 static Int ratio(const Factor& f,Int n,Int a,Int b){assert(n<Int(f.fact.size())||Int(f.fact.size())==f.modulus);Int exponent=0,out=1;bool negative=false;while(n){out=out*f.fact[n%f.modulus]%f.modulus;out=out*f.inv[a%f.modulus]%f.modulus;out=out*f.inv[b%f.modulus]%f.modulus;if((n/f.modulus-a/f.modulus-b/f.modulus)&1)negative=!negative;n/=f.p;a/=f.p;b/=f.p;exponent+=n-a-b;if(exponent>=f.q)return 0;}if(negative&&!(f.p==2&&f.q>=3))out=f.modulus-out;return out*power(f.p,exponent,f.modulus)%f.modulus;}
public:
 CombinationArbitraryMod(Int max_N,Int mod):modulus(mod){assert(max_N>=0&&1<=mod&&mod<(Int(1)<<30));Int remaining=mod;for(Int p=2;remaining>1;++p){if(p>remaining/p)p=remaining;if(remaining%p)continue;Factor f{};f.p=p;while(remaining%p==0){remaining/=p;f.modulus*=p;++f.q;}Int limit=std::min(max_N,f.modulus-1);f.fact.resize(limit+1);f.inv.resize(limit+1);f.fact[0]=1;for(Int i=1;i<=limit;++i)f.fact[i]=i%p?Int(f.fact[i-1])*i%f.modulus:f.fact[i-1];Int phi=f.modulus/p*(p-1);f.inv[limit]=power(f.fact[limit],phi-1,f.modulus);for(Int i=limit;i>0;--i)f.inv[i-1]=i%p?Int(f.inv[i])*i%f.modulus:f.inv[i];Int other=mod/f.modulus;f.coefficient=other*power(other%f.modulus,phi-1,f.modulus);factors.push_back(std::move(f));}}
 Int ncr(Int n,Int r)const{if(n<0||r<0||n<r)return 0;Int out=0;for(auto& f:factors)out=(out+ratio(f,n,r,n-r)*f.coefficient)%modulus;return out;}
 Int npr(Int n,Int r)const{if(n<0||r<0||n<r)return 0;Int out=0;for(auto& f:factors)out=(out+ratio(f,n,n-r,0)*f.coefficient)%modulus;return out;}
 Int nhr(Int n,Int r)const{if(n<0||r<0)return 0;if(!r)return 1%modulus;if(!n)return 0;assert(n<=std::numeric_limits<Int>::max()-(r-1));return ncr(n+r-1,r);}
};
inline auto initCombinationArbitraryMod(Int n,Int modulus){return CombinationArbitraryMod(n,modulus);}
inline Int ncr(const CombinationArbitraryMod& c,Int n,Int r){return c.ncr(n,r);}inline Int npr(const CombinationArbitraryMod& c,Int n,Int r){return c.npr(n,r);}inline Int nhr(const CombinationArbitraryMod& c,Int n,Int r){return c.nhr(n,r);}
}
