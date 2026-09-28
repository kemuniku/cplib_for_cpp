#pragma once
#include <cplib/convolution/convolution.hpp>
#include <optional>
#include <memory>
namespace cplib {
// x^nで切り詰め、必要なら零を補う。
template<class T> std::vector<T> prefix(const std::vector<T>& f,Int n){if(n<=0)return {};std::vector<T> out(n);std::copy_n(f.begin(),std::min(f.size(),std::size_t(n)),out.begin());return out;}
template<Modint T> T coefficient(const std::vector<T>& f,Int degree){return degree<0||degree>=Int(f.size())?T(0):f[degree];}
template<Modint T> std::vector<T> normalized(std::vector<T> f){while(!f.empty()&&f.back().val()==0)f.pop_back();return f;}
template<Modint T> std::vector<T> operator+(const std::vector<T>& f,const std::vector<T>& g){auto out=prefix(f,std::max(f.size(),g.size()));for(std::size_t i=0;i<g.size();++i)out[i]+=g[i];return out;}
template<Modint T,std::integral I> std::vector<T> operator+(std::vector<T> f,I c){if(f.empty())f.push_back(T(c));else f[0]+=T(c);return f;}
template<Modint T,std::integral I> auto operator+(I c,const std::vector<T>& f){return f+c;}
template<Modint T> std::vector<T> operator-(const std::vector<T>& f,const std::vector<T>& g){auto out=prefix(f,std::max(f.size(),g.size()));for(std::size_t i=0;i<g.size();++i)out[i]-=g[i];return out;}
template<Modint T> std::vector<T> operator-(std::vector<T> f){for(auto& x:f)x=-x;return f;}
template<Modint T> auto operator*(const std::vector<T>& f,const std::vector<T>& g){return convolution(f,g);}
template<Modint T> std::vector<T> operator*(std::vector<T> f,std::type_identity_t<T> c){for(auto& x:f)x*=c;return f;}
template<Modint T> auto operator*(std::type_identity_t<T> c,const std::vector<T>& f){return f*c;}
template<Modint T> auto operator/(const std::vector<T>& f,std::type_identity_t<T> c){return f*c.inv();}
template<Modint T> auto& operator+=(std::vector<T>& f,const std::vector<T>& g){if(f.size()<g.size())f.resize(g.size());for(std::size_t i=0;i<g.size();++i)f[i]+=g[i];return f;}
template<Modint T,std::integral I> auto& operator+=(std::vector<T>& f,I c){if(f.empty())f.push_back(T(c));else f[0]+=T(c);return f;}
template<Modint T> auto& operator-=(std::vector<T>& f,const std::vector<T>& g){if(f.size()<g.size())f.resize(g.size());for(std::size_t i=0;i<g.size();++i)f[i]-=g[i];return f;}
template<Modint T> auto& operator*=(std::vector<T>& f,std::type_identity_t<T> c){for(auto& x:f)x*=c;return f;}
template<Modint T> auto& operator/=(std::vector<T>& f,std::type_identity_t<T> c){return f*=c.inv();}
template<Modint T> std::vector<T> derivative(const std::vector<T>& f){if(f.size()<=1)return {};std::vector<T> out(f.size()-1);for(std::size_t i=1;i<f.size();++i)out[i-1]=f[i]*T(i);return out;}
namespace detail {
template<Modint T> std::vector<T> fpsIndexInverses(Int n){std::vector<T> out(n);if(n>1)out[1]=1;if(n<=2)return out;Int p=T::umod();if(!isprime(p)){for(Int i=2;i<n;++i)out[i]=T(i).inv();}else for(Int i=2;i<n;++i)out[i]=-out[p%i]*T(p/i);return out;}
// 固定側NTTを保持する特化版と巡回畳み込み版の両経路を維持。O(m log m)。
template<Modint T> void fpsInvExtend(const std::vector<T>& f,std::vector<T>& g,Int n){Int m=g.size(),size=m*2,count=n-m;auto p=T::umod();if(size>=128&&p<(1u<<30)&&(p-1)%std::uint32_t(size)==0&&isprime(Int(p))){std::vector<std::uint32_t> fixed(m),product(size);for(Int i=0;i<m;++i)fixed[i]=g[i].val();auto context=std::make_unique<avx2_ntt::FixedConvolution>(fixed.data(),m,size,p,0);Int length=std::min(Int(f.size()),n);std::vector<std::uint32_t> input(length);for(Int i=0;i<length;++i)input[i]=f[i].val();context->run(product.data(),input.data(),length);context->run(product.data(),product.data()+m,count);g.resize(n);for(Int i=0;i<count;++i)g[m+i]=-T(product[i]);}else{auto fp=f.size()<=std::size_t(n)?f:std::vector<T>(f.begin(),f.begin()+n);auto product=convolutionCyclicPowerOfTwo(fp,g,size);std::vector<T> error(count);for(Int i=0;i<count;++i)error[i]=-product[m+i];auto extension=g*error;g.resize(n);for(Int i=0;i<count;++i)g[m+i]=extension[i];}}
// Tonelli–Shanks。法は素数。
template<Modint T> std::optional<T> modSqrt(T a){Int p=T::umod();if(a.val()==0)return T(0);if(p==2)return a;if(a.pow((p-1)/2).val()!=1)return {};if(p%4==3)return a.pow((p+1)/4);Int q=p-1,s=0;while(!(q&1)){q>>=1;++s;}T z=2;while(z.pow((p-1)/2).val()!=p-1)z+=T(1);T c=z.pow(q),x=a.pow((q+1)/2),t=a.pow(q);Int m=s;while(t.val()!=1){Int i=1;T tt=t*t;while(i<m&&tt.val()!=1){tt*=tt;++i;}if(i==m)return {};T b=c.pow(Int(1)<<(m-i-1));x*=b;c=b*b;t*=c;m=i;}return x;}
}
// 素数法でO(n)、合成数法では添字ごとの逆元を使いO(n log mod)。
template<Modint T> std::vector<T> integral(const std::vector<T>& f){assert(f.size()<T::umod());auto inverses=detail::fpsIndexInverses<T>(f.size()+1);std::vector<T> out(f.size()+1);for(std::size_t i=0;i<f.size();++i)out[i+1]=f[i]*inverses[i+1];return out;}
template<Modint T> std::vector<T> inv(const std::vector<T>& f,Int n){if(n<=0)return {};assert(!f.empty()&&f[0].val()!=0);std::vector<T> out={f[0].inv()};while(Int(out.size())<n)detail::fpsInvExtend(f,out,std::min(Int(out.size())*2,n));return out;}
template<Modint T> auto inv(const std::vector<T>& f){return inv(f,Int(f.size()));}
template<Modint T> std::vector<T> operator/(const std::vector<T>& f,const std::vector<T>& g){Int n=std::max(f.size(),g.size());if(!n)return {};return prefix(f*inv(g,n),n);}
template<Modint T> auto& operator/=(std::vector<T>& f,const std::vector<T>& g){return f=f/g;}
template<Modint T> std::vector<T> log(const std::vector<T>& f,Int n){if(n<=0)return {};assert(n<=T::umod()&&!f.empty()&&f[0].val()==1);auto fp=f.size()<=std::size_t(n)?f:std::vector<T>(f.begin(),f.begin()+n);return integral(prefix(derivative(fp)*inv(fp,n-1),n-1));}
template<Modint T> auto log(const std::vector<T>& f){return log(f,Int(f.size()));}
// Newton法の誤差高次部分を巡回積で抽出して拡張する。O(n log n)。
template<Modint T> std::vector<T> exp(const std::vector<T>& f,Int n){if(n<=0)return {};assert(n<=T::umod()&&(f.empty()||f[0].val()==0));auto inverses=detail::fpsIndexInverses<T>(n);auto fp=f.size()<=std::size_t(n)?f:std::vector<T>(f.begin(),f.begin()+n);auto df=derivative(fp);std::vector<T> out={T(1)},inverse={T(1)};for(Int m=1;m<n;){Int next=std::min(m*2,n);if(Int(inverse.size())<m)detail::fpsInvExtend(out,inverse,m);auto dfp=df.size()<std::size_t(next)?df:std::vector<T>(df.begin(),df.begin()+next-1);auto product=convolutionCyclicPowerOfTwo(dfp,out,m*2);std::vector<T> error(next-m);for(Int i=0;i<Int(error.size());++i)error[i]=-product[m-1+i];auto quotient=error*prefix(inverse,error.size());for(Int i=0;i<Int(error.size());++i)error[i]=-quotient[i]*inverses[m+i];auto extension=out*error;out.resize(next);for(Int i=0;i<Int(error.size());++i)out[m+i]=extension[i];m=next;}return out;}
template<Modint T> auto exp(const std::vector<T>& f){return exp(f,Int(f.size()));}
template<Modint T> std::vector<T> pow(const std::vector<T>& f,Int k,Int n){assert(k>=0);if(n<=0)return {};assert(n<=T::umod());if(k==0){std::vector<T> out(n);out[0]=1;return out;}if(k==1)return prefix(f,n);Int ord=0;while(ord<Int(f.size())&&f[ord].val()==0)++ord;if(ord==Int(f.size())||ord>(n-1)/k)return std::vector<T>(n);Int shift=ord*k,size=n-shift;T c=f[ord],cinv=c.inv();std::vector<T> unit(std::min(Int(f.size())-ord,size));for(Int i=0;i<Int(unit.size());++i)unit[i]=f[ord+i]*cinv;auto body=exp(log(unit,size)*T(k),size);body*=c.pow(k);std::vector<T> out(n);std::copy(body.begin(),body.end(),out.begin()+shift);return out;}
template<Modint T> auto pow(const std::vector<T>& f,Int k){return pow(f,k,Int(f.size()));}
template<Modint T> std::optional<std::vector<T>> sqrt(const std::vector<T>& f,Int n){if(n<=0)return std::vector<T>{};Int ord=0;while(ord<std::min(Int(f.size()),n)&&f[ord].val()==0)++ord;if(ord==std::min(Int(f.size()),n))return std::vector<T>(n);if(ord&1)return {};Int shift=ord/2;auto root0=detail::modSqrt(f[ord]);if(!root0)return {};Int size=n-shift;std::vector<T> unit(f.begin()+ord,f.begin()+ord+std::min(Int(f.size())-ord,size)),root={*root0},inverse={root0->inv()};T half=size>1?T(2).inv():T(0);for(Int m=1;m<size;){Int next=std::min(m*2,size);if(Int(inverse.size())<m)detail::fpsInvExtend(root,inverse,m);auto square=root*root;std::vector<T> error(next-m);for(Int i=0;i<Int(error.size());++i)error[i]=(coefficient(unit,m+i)-coefficient(square,m+i))*half;auto extension=error*prefix(inverse,error.size());root.resize(next);for(Int i=0;i<Int(error.size());++i)root[m+i]=extension[i];m=next;}std::vector<T> out(n);std::copy(root.begin(),root.end(),out.begin()+shift);return out;}
template<Modint T> auto sqrt(const std::vector<T>& f){return sqrt(f,Int(f.size()));}
template<class T> struct PolynomialDivMod {std::vector<T> q,r;};
template<Modint T> PolynomialDivMod<T> divmod(const std::vector<T>& f,const std::vector<T>& g){auto a=normalized(f),b=normalized(g);assert(!b.empty());if(a.size()<b.size())return {{},std::move(a)};Int qlen=a.size()-b.size()+1;auto ar=a,br=b;std::reverse(ar.begin(),ar.end());std::reverse(br.begin(),br.end());auto q=prefix(ar*inv(br,qlen),qlen);std::reverse(q.begin(),q.end());q=normalized(std::move(q));auto r=normalized(prefix(a-q*b,b.size()-1));return {std::move(q),std::move(r)};}
template<Modint T> auto div(const std::vector<T>& f,const std::vector<T>& g){return divmod(f,g).q;}
template<Modint T> auto mod(const std::vector<T>& f,const std::vector<T>& g){return divmod(f,g).r;}
template<Modint T> auto operator%(const std::vector<T>& f,const std::vector<T>& g){return mod(f,g);}
template<Modint T> T eval(const std::vector<T>& f,std::type_identity_t<T> x){T out=0;for(Int i=Int(f.size())-1;i>=0;--i)out=out*x+f[i];return out;}
}
