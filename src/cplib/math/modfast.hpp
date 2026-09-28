#pragma once
#include <cplib/math/primitive_root.hpp>
#include <cplib/modint/modint.hpp>
#include <random>
#include <stdexcept>
namespace cplib {
// 素数法の逆元・離散対数・累乗を前計算後O(1)で返す。領域O(p^(2/3))。
template<Modint T> class ModFast {
    Int modulus=0,generator=0,bucketShift=0,powerShift=0;
    std::vector<std::pair<std::uint16_t,std::uint16_t>> fractions;
    std::vector<std::uint32_t> logarithms,powerLow,powerHigh;
    void checkModulus()const{assert(modulus==T::umod());}
    Int powerRaw(Int e)const{return Int(powerLow[e&((Int(1)<<powerShift)-1)])*Int(powerHigh[e>>powerShift])%modulus;}
    Int buildFractions(){Int n=1;while(n*n*n<modulus){n*=2;++bucketShift;}fractions.resize((modulus>>bucketShift)+1);Int a=0,b=1,c=1,d=n;while(c<=n){Int left=(a*modulus/b)>>bucketShift,right=(c*modulus/d)>>bucketShift;auto f=b<=d?std::pair{std::uint16_t(a),std::uint16_t(b)}:std::pair{std::uint16_t(c),std::uint16_t(d)};for(Int j=left;j<=right;++j)fractions[j]=f;if(c==d)break;Int k=(n+b)/d,na=c,nb=d,nc=k*c-a,nd=k*d-b;a=na;b=nb;c=nc;d=nd;}Int result=0;for(Int j=0;j<Int(fractions.size());++j){auto [fa,fb]=fractions[j];Int lo=std::max<Int>(1,j<<bucketShift),hi=std::min(modulus-1,((j+1)<<bucketShift)-1);if(lo>hi)continue;result=std::max(result,Int(fb));result=std::max(result,std::abs(lo*fb-modulus*fa));result=std::max(result,std::abs(hi*fb-modulus*fa));}return result;}
    void buildLogarithms(Int limit){
        Int p=modulus,order=p-1;logarithms.resize(limit+1);std::vector<std::uint32_t> least(limit+1);for(Int i=2;i<=limit;++i){if(least[i])continue;least[i]=i;if(i*i<=limit)for(Int j=i*i;j<=limit;j+=i)if(!least[j])least[j]=i;}
        Int step=std::min(order,Int(4)<<powerShift),capacity=1;while(capacity<2*step)capacity*=2;std::vector<std::uint32_t> keys(capacity),values(capacity);auto slot=[&](Int x){return Int((std::uint32_t(x)*2654435761u)&std::uint32_t(capacity-1));};Int value=1;
        for(Int e=0;e<step;++e){Int h=slot(value);while(keys[h])h=(h+1)&(capacity-1);keys[h]=value;values[h]=e;value=value*generator%p;}
        Int giant=powerRaw((order-step)%order);auto bsgs=[&](Int x)->Int{Int offset=0;while(offset<order){Int h=slot(x);while(keys[h]){if(Int(keys[h])==x)return (offset+values[h])%order;h=(h+1)&(capacity-1);}x=x*giant%p;offset+=step;}throw std::invalid_argument("invalid modulus or primitive root");};
        std::mt19937_64 rng(20260914);
        for(Int i=2;i<=limit;++i){if(i*i>p)logarithms[i]=(Int(logarithms[p%i])+order/2+order-Int(logarithms[p/i]))%order;else if(Int(least[i])<i)logarithms[i]=(Int(logarithms[least[i]])+Int(logarithms[i/least[i]]))%order;else if(i<100)logarithms[i]=bsgs(i);else{bool found=false;for(int attempt=0;attempt<128;++attempt){Int e=std::uniform_int_distribution<Int>(0,order-1)(rng),x=i*powerRaw(e)%p,answer=order-e;for(Int q:{2,3,5,7,11,13,17,19})while(x%q==0){x/=q;answer+=logarithms[q];}if(x>limit)continue;while(x>=i&&Int(least[x])<i){Int q=least[x];x/=q;answer+=logarithms[q];}if(x<i){answer+=logarithms[x];logarithms[i]=answer%order;found=true;break;}}if(!found)logarithms[i]=bsgs(i);}}
    }
public:
    // 構築はFarey表・篩と小素数の対数。乱択128回で打ち切り、BSGSへ移る。
    ModFast():modulus(T::umod()){assert(2<=modulus&&modulus<(Int(1)<<30));generator=modulus==2?1:primitive_root(modulus);Int width=1;while(width*width<modulus-1){width*=2;++powerShift;}powerLow.resize(width);powerHigh.resize((modulus-2)/width+1);Int value=1;for(Int i=0;i<width;++i){powerLow[i]=value;value=value*generator%modulus;}Int stride=value;value=1;for(auto& x:powerHigh){x=value;value=value*stride%modulus;}if(modulus<=64){fractions.resize(modulus);logarithms.resize(modulus);value=1;for(Int e=0;e<modulus-1;++e){fractions[value]={0,1};logarithms[value]=e;value=value*generator%modulus;}}else buildLogarithms(buildFractions());}
    Int p()const{return modulus;}
    T root()const{checkModulus();return T(generator);}
    T powRoot(Int exponent)const{checkModulus();Int e=exponent%(modulus-1);if(e<0)e+=modulus-1;return T(powerRaw(e));}
    Int log(Int x)const{checkModulus();assert(1<=x&&x<modulus);auto f=fractions[x>>bucketShift];Int a=x*Int(f.second)-modulus*Int(f.first),result=Int(logarithms[std::abs(a)])-Int(logarithms[f.second]);if(a<0)result+=(modulus-1)/2;if(result<0)result+=modulus-1;if(result>=modulus-1)result-=modulus-1;return result;}
    Int log(T x)const{return log(x.val());}
    T inv(Int x)const{Int e=log(x);return T(powerRaw(e==0?0:modulus-1-e));}
    T inv(T x)const{return inv(x.val());}
    T pow(Int x,Int exponent)const{checkModulus();assert(0<=x&&x<modulus);if(x==0){assert(exponent>=0);return T(Int(exponent==0));}Int e=exponent%(modulus-1);if(e<0)e+=modulus-1;return T(powerRaw(log(x)*e%(modulus-1)));}
    T pow(T x,Int exponent)const{return pow(x.val(),exponent);}
};
template<Modint T> auto initModFast(){return ModFast<T>();}
template<Modint T> Int p(const ModFast<T>& t){return t.p();}
template<Modint T> T root(const ModFast<T>& t){return t.root();}
template<Modint T> T powRoot(const ModFast<T>& t,Int e){return t.powRoot(e);}
template<Modint T,class X> Int log(const ModFast<T>& t,X x){return t.log(x);}
template<Modint T,class X> T inv(const ModFast<T>& t,X x){return t.inv(x);}
template<Modint T,class X> T pow(const ModFast<T>& t,X x,Int e){return t.pow(x,e);}
}
