#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 各接尾辞と全体の共通接頭辞長を求める。O(n)。
template<class Range> std::vector<Int> zalgorithm(const Range& s){Int n=s.size();std::vector<Int> result(n,-1);if(!n)return result;result[0]=n;Int i=1,j=0;while(i<n){while(i+j<n&&s[j]==s[i+j])++j;result[i]=j;if(!j){++i;continue;}Int k=1;while(i+k<n&&k+result[k]<j){result[i+k]=result[k];++k;}i+=k;j-=k;}return result;}
}
