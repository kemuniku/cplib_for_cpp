#pragma once
#include <cplib/modint/modint.hpp>
namespace cplib {
// 最小線形漸化式。末尾の零も次数の一部として残す。O(N²)、領域O(N)。
template<Modint T> std::vector<T> berlekampMassey(const std::vector<T>& a){std::vector<T> connection={T(1)},previous={T(1)};Int order=0,shift=1;T previousDiscrepancy=1;for(Int n=0;n<Int(a.size());++n){T discrepancy=a[n];for(Int i=1;i<=order;++i)discrepancy+=connection[i]*a[n-i];if(discrepancy.val()==0){++shift;continue;}auto old=connection;T scale=discrepancy/previousDiscrepancy;if(connection.size()<previous.size()+shift)connection.resize(previous.size()+shift);for(Int i=0;i<Int(previous.size());++i)connection[i+shift]-=scale*previous[i];if(2*order<=n){order=n+1-order;previous=std::move(old);previousDiscrepancy=discrepancy;shift=1;}else ++shift;}std::vector<T> out(order);for(Int i=0;i<order;++i)out[i]=-connection[i+1];return out;}
}
