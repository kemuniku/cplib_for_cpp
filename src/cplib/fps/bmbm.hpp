#pragma once
#include <cplib/fps/berlekamp_massey.hpp>
#include <cplib/fps/bostan_mori.hpp>
namespace cplib {
// 漸化式を推定して第k項を計算。O(N²+M(D)log(k+1))。
template<Modint T> T bmbm(const std::vector<T>& a,Int k){assert(k>=0);if(k<Int(a.size()))return a[k];auto coefficients=berlekampMassey(a);if(coefficients.empty())return T(0);return linearRecurrenceKth(std::vector<T>(a.begin(),a.begin()+coefficients.size()),coefficients,k);}
}
