#ifndef CPLIB_HASH_HEADER
#define CPLIB_HASH_HEADER <cplib/str/hash_string.hpp>
#endif
#include CPLIB_HASH_HEADER
#include <random>
using namespace cplib;
int main(){std::mt19937 rng(452);assert(tohash(65)==tohash('A'));
    for(Int n:{0,1,2,5,64,500,501001}){std::string text(n,'a');for(char& c:text)c=char(rng()%256);auto h=initRollingHash(text);auto rev=initRollingHash(text);assert(toHashString(h)==tohash(text)&&toHashString(rev)==cplib::tohash(text));for(int q=0;q<200;++q){Int l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);auto s=h[closed_slice(l,r-1)];auto t=rev[closed_slice(l,r-1)];auto sub=text.substr(l,r-l);assert(s.to_string()==sub&&toHashString(s)==tohash(sub)&&toHashString(t)==cplib::tohash(sub));auto reverse=sub;std::reverse(reverse.begin(),reverse.end());
#ifdef CPLIB_REVERSIBLE_HASH
assert(reversed(toHashString(t))==cplib::tohash(reverse));
#endif

#ifdef CPLIB_REVERSIBLE_HASH
assert(t.isPalindrome()==(sub==reverse));
#endif
Int split=rng()%(sub.size()+1);auto a=s.substr(0,split),b=s.substr(split,s.len());assert((a&b)==toHashString(s));
#ifndef CPLIB_REVERSIBLE_HASH
assert(removePrefix(toHashString(s),toHashString(a))==toHashString(b));
#endif
auto ra=t.substr(0,split),rb=t.substr(split,t.len());assert((ra&rb)==toHashString(t));if(q<10&&sub.size()<1000){Int k=rng()%10;std::string repeated;for(Int i=0;i<k;++i)repeated+=sub;assert(toHashString(s)*k==tohash(repeated));auto value=toHashString(t)*k;assert(value==cplib::tohash(repeated));auto reversed_text=repeated;std::reverse(reversed_text.begin(),reversed_text.end());
#ifdef CPLIB_REVERSIBLE_HASH
assert(reversed(value)==cplib::tohash(reversed_text));
#endif
}}}
    for(int q=0;q<1000;++q){std::string a(rng()%100,'a'),b(rng()%100,'a');for(char& c:a)c='a'+rng()%3;for(char& c:b)c='a'+rng()%3;auto x=initRollingHash(a),y=initRollingHash(b);auto rx=cplib::initRollingHash(a),ry=cplib::initRollingHash(b);Int common=0;while(common<Int(std::min(a.size(),b.size()))&&a[common]==b[common])++common;assert(LCP(x,y)==common&&LCP(rx,ry)==common&&(x<y)==(a<b)&&(rx<ry)==(a<b));}
    auto hash=tohash(std::string("abc"));assert((hash*1000000000000LL).len()==3000000000000LL);assert(hash*100==((hash*50)&(hash*50)));
}
