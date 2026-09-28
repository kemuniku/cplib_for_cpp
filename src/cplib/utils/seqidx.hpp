#pragma once
#include <cplib/common.hpp>
#include <functional>
#include <iterator>
namespace cplib {
// Nimのidx,itを引数(idx,it)で受け取る。いずれもO(N)。
template<class Range,class Pred> bool allItidx(const Range& s,Pred pred){Int i=0;for(const auto& x:s)if(!pred(i++,x))return false;return true;}
template<class Range,class Pred> Int findItIdx(const Range& s,Pred pred){Int i=0;for(const auto& x:s){if(pred(i,x))return i;++i;}return -1;}
template<class Range,class Pred> bool anyItIdx(const Range& s,Pred pred){return findItIdx(s,pred)!=-1;}
template<class Range,class Op> auto mapItIdx(const Range& s,Op op){using T=std::decay_t<std::invoke_result_t<Op,Int,decltype(*std::begin(s))>>;std::vector<T> result;result.reserve(std::size(s));Int i=0;for(const auto& x:s)result.push_back(op(i++,x));return result;}
template<class Init> auto newSeqWithIdx(Int n,Init init){assert(n>=0);using T=std::decay_t<std::invoke_result_t<Init,Int>>;std::vector<T> result;result.reserve(n);for(Int i=0;i<n;++i)result.push_back(init(i));return result;}
}
