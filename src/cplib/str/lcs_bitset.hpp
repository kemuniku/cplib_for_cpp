#pragma once
#include <cplib/collections/bitset_avx512.hpp>
#include <unordered_map>
namespace cplib {
namespace detail {
// vector<bool>や独自の等値比較要素もコピーせず区間を参照する。
template<class Sequence> struct View {
 using value_type=typename Sequence::value_type;
 const Sequence* source;std::size_t begin,end;
 std::size_t size()const{return end-begin;}
 decltype(auto) operator[](std::size_t i)const{return (*source)[begin+i];}
 View slice(std::size_t l,std::size_t r)const{return {source,begin+l,begin+r};}
};
template<bool Backwards=false,class A,class B> BitSetAvx512 dp(const A& a,const B& b){
 using T=typename B::value_type;
 Int m=b.size(),words=(m+63)/64;BitSetAvx512 state(m),matched(m);
 auto advance=[&](const BitSetAvx512& mask){cplib::fuse([&](auto& f){auto s=f.var(state);auto equal=f.read(mask);auto combined=s|equal;auto difference=combined+~(s<<1);s=combined&~difference;});};
 auto av=[&](std::size_t i)->decltype(auto){return a[Backwards?a.size()-1-i:i];};
 auto bv=[&](std::size_t i)->decltype(auto){return b[Backwards?b.size()-1-i:i];};
 if constexpr(std::is_integral_v<T>||std::is_enum_v<T>||std::is_same_v<T,std::string>){
  std::unordered_map<T,Int> ids;std::vector<std::vector<Int>> positions;
  for(Int i=0;i<m;++i){auto [it,inserted]=ids.try_emplace(bv(i),positions.size());if(inserted)positions.emplace_back();positions[it->second].push_back(i);}
  std::vector<BitSetAvx512> masks(positions.size());
  // 頻度が語数を超える種類だけを事前にビット列化し、追加領域O(m)を保つ。
  for(std::size_t id=0;id<positions.size();++id)if(Int(positions[id].size())>words)masks[id]=cplib::initBitSetFromIndexes(positions[id],m);
  for(std::size_t k=0;k<a.size();++k){auto it=ids.find(av(k));if(it==ids.end())continue;Int id=it->second;if(masks[id].len())advance(masks[id]);else{matched.clear();for(Int i:positions[id])matched[i]=true;advance(matched);}}
 }else{
  for(std::size_t k=0;k<a.size();++k){matched.clear();for(Int i=0;i<m;++i)if(bv(i)==av(k))matched[i]=true;advance(matched);}
 }
 return state;
}
template<class A,class B> std::size_t split(const A& a,const B& b,std::size_t middle){
 auto forward=dp(a.slice(0,middle),b),backward=dp<true>(a.slice(middle,a.size()),b);
 Int left=0,right=backward.popcount(),best=right;std::size_t result=0;
 for(std::size_t j=0;j<b.size();++j){left+=bool(forward[j]);right-=bool(backward[b.size()-1-j]);if(left+right>best){best=left+right;result=j+1;}}
 return result;
}
template<class A,class B,class T> void restore(const A& a,const B& b,std::vector<T>& output){
 if(a.size()<b.size()){restore(b,a,output);return;}if(!b.size())return;
 if(b.size()==1){for(std::size_t i=0;i<a.size();++i)if(a[i]==b[0]){output.push_back(b[0]);break;}return;}
 auto middle=a.size()/2;auto cut=split(a,b,middle);
 if(cut)restore(a.slice(0,middle),b.slice(0,cut),output);
 if(cut<b.size())restore(a.slice(middle,a.size()),b.slice(cut,b.size()),output);
}
}
// 整数・文字・bool・enum・stringは期待O(m+n ceil(m/64))、追加領域O(m)。
// その他は等値比較だけを要求しO(nm)。n>=m、空入力はO(1)。
template<class A,class B> Int LCS(const A& a,const B& b){if(a.size()<b.size())return cplib::LCS(b,a);if(!b.size())return 0;return detail::dp(detail::View<A>{&a,0,a.size()},detail::View<B>{&b,0,b.size()}).popcount();}
// Hirschberg法。ハッシュ可能な基本型は期待O(nm/64+(n+m)log(n+m))、領域O(n+m)。
template<class A,class B> auto restoreLCS(const A& a,const B& b){std::vector<typename A::value_type> output;detail::restore(detail::View<A>{&a,0,a.size()},detail::View<B>{&b,0,b.size()},output);return output;}
}
