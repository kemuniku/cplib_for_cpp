#pragma once
#include <cplib/math/isqrt.hpp>
#include <cplib/math/isprime.hpp>
#include <cplib/modint/modint.hpp>
namespace cplib {
namespace detail {
template<class T> T multiplicative_polynomial(const std::vector<T>& c,Int n){T out=0,x=n;for(auto it=c.rbegin();it!=c.rend();++it)out=out*x+*it;return out;}
template<class T> auto multiplicative_power_sums(Int degree){std::vector<std::vector<T>> out;for(Int j=0;j<=degree;++j){std::vector<T> binomial(j+2);binomial[0]=1;for(Int k=1;k<=j+1;++k)binomial[k]=binomial[k-1]*T(j+2-k)/T(k);auto poly=binomial;poly[0]-=1;for(Int k=0;k<j;++k)for(Int t=0;t<Int(out[k].size());++t)poly[t]-=binomial[k]*out[k][t];T inv=T(j+1).inv();for(auto& x:poly)x*=inv;out.push_back(std::move(poly));}return out;}
}
template<class T,class PrimePower> T multiplicativePrefixSum(Int n,std::span<const T> primeCoefficients,PrimePower primePower){
 assert(n>=0&&n<=(Int(1)<<40)&&isprime(T::umod()));if(n<=1)return T(n);std::vector<T> c(primeCoefficients.begin(),primeCoefficients.end());while(!c.empty()&&c.back().val()==0)c.pop_back();
 if(n<=4096){std::vector<Int> least(n+1);std::vector<T> values(n+1);values[1]=1;T out=1;for(Int p=2;p<=n;++p){if(!least[p])for(Int m=p;m<=n;m+=p)if(!least[m])least[m]=p;Int rest=p,prime=least[p],exponent=0;while(rest%prime==0){rest/=prime;++exponent;}values[p]=values[rest]*(exponent==1?detail::multiplicative_polynomial(c,prime):primePower(prime,exponent));out+=values[p];}return out;}
 assert(c.size()<T::umod());Int root=isqrt(n),largeCount=n/(root+1);std::vector<bool> composite(root+1);std::vector<Int> primes;for(Int p=2;p<=root;++p){if(composite[p])continue;primes.push_back(p);if(p<=root/p)for(Int m=p*p;m<=root;m+=p)composite[m]=true;}
 std::vector<T> small(root+1),large(largeCount+1);auto powerSums=detail::multiplicative_power_sums<T>(Int(c.size())-1);
 for(Int d=0;d<Int(c.size());++d){if(c[d].val()==0)continue;std::vector<T> low(root+1),high(largeCount+1);for(Int x=1;x<=root;++x)low[x]=detail::multiplicative_polynomial(powerSums[d],x)-1;for(Int i=1;i<=largeCount;++i)high[i]=detail::multiplicative_polynomial(powerSums[d],n/i)-1;for(Int p:primes){T weight=T(p).pow(d),before=low[p-1];Int bound=std::min(largeCount,n/(p*p)),ip=p;for(Int i=1;i<=bound;++i,ip+=p)high[i]-=weight*((ip<=largeCount?high[ip]:low[n/ip])-before);for(Int x=root;x>=p*p;--x)low[x]-=weight*(low[x/p]-before);}for(Int x=1;x<=root;++x)small[x]+=c[d]*low[x];for(Int i=1;i<=largeCount;++i)large[i]+=c[d]*high[i];}
 std::vector<Int> offsets(primes.size());std::vector<T> values,primePrefix(primes.size()+1);for(Int i=0;i<Int(primes.size());++i){Int p=primes[i];offsets[i]=values.size();T value=detail::multiplicative_polynomial(c,p);primePrefix[i+1]=primePrefix[i]+value;values.push_back(value);Int power=p,e=1;while(power<=n/p){power*=p;values.push_back(primePower(p,++e));}}
 T answer=T(1)+large[1];auto visit=[&](auto&& self,Int limit,Int start,T weight)->void{if(weight.val()==0)return;for(Int i=start;i<Int(primes.size())&&primes[i]*primes[i]<=limit;++i){Int p=primes[i],power=p,e=1;while(power<=limit/p){Int remaining=limit/power;T value=weight*values[offsets[i]+e-1],sum=remaining<=root?small[remaining]:large[n/remaining];answer+=weight*values[offsets[i]+e]+value*(sum-primePrefix[i+1]);self(self,remaining,i+1,value);power*=p;++e;}}};visit(visit,n,0,T(1));return answer;
}
template<class T,class F> T multiplicativePrefixSum(Int n,const std::vector<T>& c,F f){return multiplicativePrefixSum<T>(n,std::span<const T>(c),f);}
}
