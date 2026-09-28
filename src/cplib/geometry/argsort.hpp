#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 128ビット積で整数偏角を比較。正のx軸から反時計回り。
inline Int argcmp(std::pair<Int,Int> l,std::pair<Int,Int> r){bool a=l.second>0||(l.second==0&&l.first>0),b=r.second>0||(r.second==0&&r.first>0);if(a!=b)return b<a?-1:1;__int128_t x=__int128_t(l.second)*r.first,y=__int128_t(l.first)*r.second;return x<y?-1:x>y?1:0;}
inline void argsort(std::vector<std::pair<Int,Int>>& x){std::stable_sort(x.begin(),x.end(),[](auto l,auto r){return argcmp(l,r)<0;});}
inline auto argsorted(std::vector<std::pair<Int,Int>> x){argsort(x);return x;}
}
