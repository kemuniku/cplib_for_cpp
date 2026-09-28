#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <memory>
#include <random>
#include <string_view>
namespace cplib::detail {
inline constexpr UInt hash_string_mod=(UInt(1)<<61)-1;
inline UInt hash_string_modulo(UInt x){x=(x>>61)+(x&hash_string_mod);if(x>=hash_string_mod)x-=hash_string_mod;return x;}
// 31/30 bit 分割による元の乗算。128bit乗算への置き換えはしない。
inline UInt hash_string_mul(UInt a,UInt b){UInt au=a>>31,al=a&((UInt(1)<<31)-1),bu=b>>31,bl=b&((UInt(1)<<31)-1),mid=al*bu+au*bl;return au*bu*2+(mid>>30)+((mid&((UInt(1)<<30)-1))<<31)+al*bl;}
inline UInt hash_string_pow(UInt a,Int n){UInt result=1;while(n>0){if(n&1)result=hash_string_modulo(hash_string_mul(result,a));a=hash_string_modulo(hash_string_mul(a,a));n>>=1;}return result;}
template<bool Reverse> struct HashStringConfig {UInt base,inverse;std::vector<UInt> pows,invpows;HashStringConfig():pows(500001),invpows(500001){std::mt19937_64 rng(std::random_device{}());base=std::uniform_int_distribution<UInt>(129,UInt(1)<<30)(rng);inverse=hash_string_pow(base,hash_string_mod-2);pows[0]=invpows[0]=1;for(std::size_t i=1;i<pows.size();++i){pows[i]=hash_string_modulo(hash_string_mul(pows[i-1],base));invpows[i]=hash_string_modulo(hash_string_mul(invpows[i-1],inverse));}}
    UInt power(Int n)const{assert(n>=0);return n<Int(pows.size())?pows[n]:hash_string_pow(base,n);}UInt invpower(Int n)const{assert(n>=0);return n<Int(invpows.size())?invpows[n]:hash_string_pow(inverse,n);}};
template<bool R> inline HashStringConfig<R> hash_string_config;
template<bool Reverse> struct ReverseHashField {};template<> struct ReverseHashField<true>{UInt rhash=0;};
}
namespace cplib {
template<bool Reverse> struct BasicHashString:detail::ReverseHashField<Reverse> {
    UInt hash=0,bpow=1;Int size=0;
    static BasicHashString character(Int s){if constexpr(Reverse)assert(0<=s&&UInt(s)<detail::hash_string_mod);BasicHashString out;out.hash=UInt(s)%detail::hash_string_mod;if constexpr(Reverse)out.rhash=out.hash;out.bpow=detail::hash_string_config<Reverse>.base;out.size=1;return out;}
    template<class Range> static BasicHashString sequence(const Range& s){BasicHashString out;UInt tmp=1;for(Int i=Int(s.size());i-->0;){auto value=[](auto x){if constexpr(std::is_same_v<decltype(x),char>)return Int(static_cast<unsigned char>(x));else return Int(x);};out.hash=detail::hash_string_modulo(out.hash+detail::hash_string_mul(character(value(s[i])).hash,tmp));if constexpr(Reverse)out.rhash=detail::hash_string_modulo(out.rhash+detail::hash_string_mul(character(value(s[s.size()-1-i])).hash,tmp));tmp=detail::hash_string_modulo(detail::hash_string_mul(tmp,detail::hash_string_config<Reverse>.base));}out.bpow=detail::hash_string_config<Reverse>.power(s.size());out.size=s.size();return out;}
    Int len()const{return size;}
    friend BasicHashString operator&(const BasicHashString& l,const BasicHashString& r){BasicHashString out;out.hash=detail::hash_string_modulo(detail::hash_string_modulo(detail::hash_string_mul(l.hash,r.bpow))+r.hash);if constexpr(Reverse)out.rhash=detail::hash_string_modulo(detail::hash_string_modulo(detail::hash_string_mul(r.rhash,l.bpow))+l.rhash);out.bpow=detail::hash_string_modulo(detail::hash_string_mul(l.bpow,r.bpow));out.size=l.size+r.size;return out;}
    bool operator==(const BasicHashString& r)const{return size==r.size&&hash==r.hash;}
    BasicHashString operator*(Int x)const{BasicHashString out,tmp=*this;Int total=size*x;while(x>0){if(x&1)out=out&tmp;if(x>1)tmp=tmp&tmp;x>>=1;}out.size=total;return out;}
    BasicHashString removePrefix(const BasicHashString& prefix)const requires (!Reverse){Int n=size-prefix.size;BasicHashString out;out.hash=detail::hash_string_modulo(hash+detail::hash_string_mod-detail::hash_string_modulo(detail::hash_string_mul(prefix.hash,detail::hash_string_config<Reverse>.power(n))));out.bpow=detail::hash_string_config<Reverse>.power(n);out.size=n;return out;}
    bool isPalindrome()const requires Reverse{return hash==this->rhash;}BasicHashString reversed()const requires Reverse{auto out=*this;std::swap(out.hash,out.rhash);return out;}
};
template<bool Reverse> struct BasicRollingHash {
    using Hash=BasicHashString<Reverse>;
    struct Base {std::string S;std::vector<UInt> prefixs,rprefixs;Int len()const{return S.size();}};
    std::shared_ptr<Base> R;Int l=0,r=0;
    BasicRollingHash()=default;explicit BasicRollingHash(std::string_view s):R(std::make_shared<Base>()),r(s.size()){R->S=s;R->prefixs.resize(s.size()+1);if constexpr(Reverse)R->rprefixs.resize(s.size()+1);UInt tmp=1;for(std::size_t i=1;i<=s.size();++i){UInt ch=static_cast<unsigned char>(s[i-1]);R->prefixs[i]=detail::hash_string_modulo(detail::hash_string_mul(R->prefixs[i-1],detail::hash_string_config<Reverse>.base)+ch);if constexpr(Reverse){R->rprefixs[i]=detail::hash_string_modulo(detail::hash_string_mul(ch,tmp)+R->rprefixs[i-1]);tmp=detail::hash_string_modulo(detail::hash_string_mul(tmp,detail::hash_string_config<Reverse>.base));}}}
    Int len()const{return r-l;}char operator[](Int i)const{assert(0<=i&&i<len());return R->S[l+i];}char operator[](BackwardsIndex i)const{return (*this)[len()-i.value];}
    BasicRollingHash substr(Int a,Int b)const{assert(0<=a&&a<=b&&b<=len());auto out=*this;if(a==b)out.l=out.r=0;else{out.l=l+a;out.r=l+b;}return out;}
    template<class A,class B> auto operator[](ClosedSlice<A,B> s)const{return substr(resolve_index(len(),s.a),resolve_index(len(),s.b)+1);}
    UInt gethash(Int a,Int b)const{return detail::hash_string_modulo(R->prefixs[l+b]+detail::hash_string_mod-detail::hash_string_modulo(detail::hash_string_mul(R->prefixs[l+a],detail::hash_string_config<Reverse>.power(b-a))));}
    operator Hash()const{Hash out;out.hash=gethash(0,len());out.bpow=detail::hash_string_config<Reverse>.power(len());out.size=len();if constexpr(Reverse)out.rhash=detail::hash_string_modulo(detail::hash_string_mul(detail::hash_string_modulo(R->rprefixs[r]+detail::hash_string_mod-R->rprefixs[l]),detail::hash_string_config<Reverse>.invpower(l)));return out;}
    bool operator==(const BasicRollingHash& t)const{return len()==t.len()&&gethash(0,len())==t.gethash(0,t.len());}
    std::string to_string()const{return R->S.substr(l,len());}
    friend Hash operator&(const BasicRollingHash& a,const BasicRollingHash& b){return Hash(a)&Hash(b);}
    bool isPalindrome()const requires Reverse{return Hash(*this).isPalindrome();}
};
template<bool R> Int len(const BasicHashString<R>& s){return s.len();}template<bool R> Int len(const BasicRollingHash<R>& s){return s.len();}
template<bool R> auto toHashString(const BasicRollingHash<R>& s){return BasicHashString<R>(s);}
template<bool R> Int LCP(const BasicRollingHash<R>& s,const BasicRollingHash<R>& t){Int ok=0,ng=std::min(s.len(),t.len())+1;while(ng-ok>1){Int mid=(ok+ng)/2;if(s.gethash(0,mid)==t.gethash(0,mid))ok=mid;else ng=mid;}return ok;}
template<bool R> int cmp(const BasicRollingHash<R>& s,const BasicRollingHash<R>& t){Int n=LCP(s,t);if(n==std::min(s.len(),t.len()))return (s.len()>t.len())-(s.len()<t.len());return static_cast<unsigned char>(s[n])<static_cast<unsigned char>(t[n])?-1:1;}
template<bool R> bool operator<(const BasicRollingHash<R>& s,const BasicRollingHash<R>& t){return cmp(s,t)<0;}
template<bool R> std::string to_string(const BasicRollingHash<R>& s){return s.to_string();}
inline auto removePrefix(const BasicHashString<false>& s,const BasicHashString<false>& p){return s.removePrefix(p);}
inline bool isPalindrome(const BasicHashString<true>& s){return s.isPalindrome();}inline auto reversed(const BasicHashString<true>& s){return s.reversed();}
inline bool isPalindrome(const BasicRollingHash<true>& s){return s.isPalindrome();}inline auto reversed(const BasicRollingHash<true>& s){return BasicHashString<true>(s).reversed();}
}
