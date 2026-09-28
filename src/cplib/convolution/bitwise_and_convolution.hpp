#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 上位集合ゼータ変換によるAND畳み込み。O(N log N)、領域O(N)。
template<class T> std::vector<T> bitwiseAndConvolution(std::vector<T> a,std::vector<T> b){assert(a.size()==b.size());std::size_t n=a.size();if(!n)return {};assert(std::has_single_bit(n));for(std::size_t bit=1;bit<n;bit<<=1)for(std::size_t mask=0;mask<n;++mask)if(!(mask&bit)){a[mask]+=a[mask|bit];b[mask]+=b[mask|bit];}for(std::size_t mask=0;mask<n;++mask)a[mask]*=b[mask];for(std::size_t bit=1;bit<n;bit<<=1)for(std::size_t mask=0;mask<n;++mask)if(!(mask&bit))a[mask]-=a[mask|bit];return a;}
}
