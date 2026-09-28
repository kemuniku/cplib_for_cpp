#pragma once
#include <cplib/str/fixedlength_merged_static_string.hpp>
namespace cplib {
template<class T> struct RepeatedStaticString {StaticString<T> period;Int size=0;Int len()const{return size;}T operator[](Int i)const{assert(0<=i&&i<size);return period[i%period.len()];}T operator[](BackwardsIndex i)const{return (*this)[size-i.value];}std::string to_string()const{std::ostringstream out;for(Int i=0;i<size;++i){if constexpr(!std::is_same_v<T,char>){if(i)out<<' ';}out<<(*this)[i];}return out.str();}};
template<class T> auto initRepeatedStaticString(const StaticString<T>& s,Int k){assert(k>=0&&(k==0||s.len()>0));return RepeatedStaticString<T>{s,k};}
namespace detail {
// 無限反復列の最初の不一致は ST と TS の最初の不一致に等しい。O(1)。
template<class T> Int infinite_lcp(const StaticString<T>& s,const StaticString<T>& t){assert(s.base==t.base&&s.len()>0&&t.len()>0);Int result=lcp(s&t,t&s);return result==s.len()+t.len()?std::numeric_limits<Int>::max():result;}
template<class S,class T> int repeated_cmp(const S& s,const T& t,Int common){if(common==std::min(s.len(),t.len()))return (s.len()>t.len())-(s.len()<t.len());return static_value_less(s[common],t[common])?-1:1;}
}
template<class T> Int lcp(const RepeatedStaticString<T>& s,const RepeatedStaticString<T>& t){assert(s.period.base==t.period.base);Int n=std::min(s.len(),t.len());return n?std::min(n,detail::infinite_lcp(s.period,t.period)):0;}
template<class T> Int lcp(const RepeatedStaticString<T>& s,const StaticString<T>& t){assert(s.period.base==t.base);Int n=std::min(s.len(),t.len());return n?std::min(n,detail::infinite_lcp(s.period,t)):0;}
template<class T> Int lcp(const StaticString<T>& s,const RepeatedStaticString<T>& t){return lcp(t,s);}
template<class T> int cmp(const RepeatedStaticString<T>& s,const RepeatedStaticString<T>& t){return detail::repeated_cmp(s,t,lcp(s,t));}template<class T> int cmp(const RepeatedStaticString<T>& s,const StaticString<T>& t){return detail::repeated_cmp(s,t,lcp(s,t));}template<class T> int cmp(const StaticString<T>& s,const RepeatedStaticString<T>& t){return detail::repeated_cmp(s,t,lcp(s,t));}
template<class T> bool operator==(const RepeatedStaticString<T>& s,const RepeatedStaticString<T>& t){return s.len()==t.len()&&lcp(s,t)==s.len();}template<class T> auto operator<=>(const RepeatedStaticString<T>& s,const RepeatedStaticString<T>& t){return cmp(s,t)<=>0;}
template<class T> bool operator==(const RepeatedStaticString<T>& s,const StaticString<T>& t){return s.len()==t.len()&&lcp(s,t)==s.len();}template<class T> auto operator<=>(const RepeatedStaticString<T>& s,const StaticString<T>& t){return cmp(s,t)<=>0;}
template<class T> Int len(const RepeatedStaticString<T>& s){return s.len();}template<class T> std::string to_string(const RepeatedStaticString<T>& s){return s.to_string();}
}
