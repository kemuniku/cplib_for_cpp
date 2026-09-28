#pragma once
#include <cplib/str/static_string.hpp>
#include <stdexcept>
namespace cplib::detail {
template<class S> Int merged_length(const S& s){Int out=0;for(std::size_t i=0;i<s.L.size();++i)out+=s.R[i]-s.L[i];return out;}
template<class S> auto merged_at(const S& s,Int idx){assert(idx>=0);for(std::size_t i=0;i<s.L.size();++i){Int length=s.R[i]-s.L[i];if(idx<length)return s.base->S[s.L[i]+idx];idx-=length;}throw std::out_of_range("merged string index");}
// 区間境界だけを進め、各区間内の LCP は RMQ で取得。O(結合数の和)。
template<class S,class T> Int merged_lcp(const S& s,const T& t){if(s.L.empty()||t.L.empty())return 0;assert(s.base==t.base);std::size_t si=0,ti=0;Int sl=s.L[0],tl=t.L[0],result=0;while(si<s.L.size()&&ti<t.L.size()){Int sn=s.R[si]-sl,tn=t.R[ti]-tl,n=static_lcp_range(s.base,sl,s.R[si],tl,t.R[ti]);result+=n;if(n<sn&&n<tn)return result;if(n==sn){if(++si==s.L.size())return result;sl=s.L[si];}else sl+=n;if(n==tn){if(++ti==t.L.size())return result;tl=t.L[ti];}else tl+=n;}return result;}
template<class S,class T> int merged_cmp(const S& s,const T& t){std::size_t si=0,ti=0;Int sl=s.L.empty()?0:s.L[0],tl=t.L.empty()?0:t.L[0];for(;;){while(si<s.L.size()&&sl==s.R[si])if(++si<s.L.size())sl=s.L[si];while(ti<t.L.size()&&tl==t.R[ti])if(++ti<t.L.size())tl=t.L[ti];if(si==s.L.size())return ti==t.L.size()?0:-1;if(ti==t.L.size())return 1;assert(s.base==t.base);Int limit=std::min(s.R[si]-sl,t.R[ti]-tl),n=static_lcp_range(s.base,sl,s.R[si],tl,t.R[ti]);if(n<limit)return static_value_less(s.base->S[sl+n],t.base->S[tl+n])?-1:1;sl+=n;tl+=n;}}
template<class S> std::string merged_string(const S& s){std::ostringstream out;bool first=true;for(std::size_t i=0;i<s.L.size();++i)for(Int j=s.L[i];j<s.R[i];++j){if constexpr(!std::is_same_v<typename S::value_type,char>){if(!first)out<<' ';}first=false;out<<s.base->S[j];}return out.str();}
}
