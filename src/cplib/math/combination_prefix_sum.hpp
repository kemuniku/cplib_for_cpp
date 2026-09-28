#pragma once
#include <cplib/math/combination.hpp>
#include <cplib/utils/mo.hpp>
#include <tuple>
namespace cplib {
// 二項係数のprefix sumをMoでオフライン処理。N=max n、Q=クエリ数としてO(N sqrt Q+Q log Q)。
template<class M> std::vector<M> combinationPrefixSum(std::span<const std::pair<Int,Int>> queries){if(queries.empty())return {};Int maxX=0;for(auto [x,y]:queries)maxX=std::max(maxX,x);auto c=initCombination<M>(std::max<Int>(1,maxX));auto mo=initMo(std::max<Int>(1,maxX),queries.size());for(auto [x,y]:queries)mo.insert(x,std::min(x,y));std::vector<M> answers(queries.size());Int x=0,y=0;M sum=1;mo.run([&](Int newX){x=newX;sum=(sum+c.ncr(x,y))/2;},[&](Int oldY){y=oldY+1;sum+=c.ncr(x,y);},[&](Int oldX){sum=2*sum-c.ncr(oldX,y);x=oldX+1;},[&](Int newY){sum-=c.ncr(x,newY+1);y=newY;},[&](Int i){answers[i]=sum;});return answers;}
template<class M> std::vector<M> combinationRangeSum(std::span<const std::tuple<Int,Int,Int>> queries){std::vector<M> result(queries.size());std::vector<std::pair<Int,Int>> prefixes,indices;for(Int i=0;i<Int(queries.size());++i){auto [n,l,r]=queries[i];Int left=std::max<Int>(0,l),right=std::min(n+1,r);if(left>=right)continue;indices.emplace_back(i,prefixes.size());prefixes.emplace_back(n,right-1);prefixes.emplace_back(n,left-1);}auto sums=combinationPrefixSum<M>(prefixes);for(auto [i,p]:indices)result[i]=sums[p]-sums[p+1];return result;}
}
