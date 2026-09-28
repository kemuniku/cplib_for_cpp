#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 奇数長回文の半径（中心込み）を求める。O(n)。
template<class Range> std::vector<Int> manacher(const Range& s){Int n=s.size(),c=0;std::vector<Int> result(n);for(Int i=0;i<n;++i){Int l=2*c-i;if(0<=l&&l<n&&i+result[l]<c+result[c])result[i]=result[l];else{Int j=c+result[c]-i;while(i-j>=0&&i+j<n&&s[i-j]==s[i+j])++j;result[i]=j;c=i;}}return result;}
// 2n-1個の中心の最大回文区間を返す。partitionは入力に含まない。O(n)。
template<class Range,class T> std::vector<std::pair<Int,Int>> get_palindromes(const Range& s,T partition){Int n=s.size();if(!n)return {};std::vector<T> tmp(2*n-1,partition);for(Int i=0;i<n;++i)tmp[2*i]=s[i];auto mana=manacher(tmp);std::vector<std::pair<Int,Int>> result(2*n-1);for(Int i=0;i<2*n-1;++i){if(i%2==0){Int x=(mana[i]+1)/2-1,idx=i/2;result[i]={idx-x,idx+x+1};}else{Int x=mana[i]/2,idx=i/2+1;result[i]=x?std::pair{idx-x,idx+x}:std::pair<Int,Int>{-1,-1};}}return result;}
inline auto get_palindromes(const std::string& s,char partition='$'){return get_palindromes<std::string,char>(s,partition);}
}
