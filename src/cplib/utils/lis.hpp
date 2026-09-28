#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 狭義単調増加部分列の最大長を返す。O(n log n)。
template<class Range> Int lis(const Range& a){
    std::vector<typename Range::value_type> dp;
    for(const auto& x:a){auto it=std::lower_bound(dp.begin(),dp.end(),x);if(it==dp.end())dp.push_back(x);else *it=x;}return dp.size();
}
// 元実装と同じ後方走査でLISの添字列を復元する。O(n log n)。
template<class Range> std::vector<Int> restore_lis_index(const Range& a){
    std::vector<Int> p(a.size()),result;std::vector<typename Range::value_type> dp;
    for(std::size_t i=0;i<a.size();++i){auto pos=std::lower_bound(dp.begin(),dp.end(),a[i])-dp.begin();if(std::size_t(pos)==dp.size())dp.push_back(a[i]);else dp[pos]=a[i];p[i]=pos;}
    Int t=Int(dp.size())-1;
    for(Int i=Int(a.size())-1;i>=0;--i)if(p[i]==t){result.push_back(i);--t;}
    std::reverse(result.begin(),result.end());return result;
}
// LISの値列を復元する。O(n log n)。
template<class Range> auto restore_lis(const Range& a){std::vector<typename Range::value_type> result;for(Int i:restore_lis_index(a))result.push_back(a[i]);return result;}
}
