#pragma once
#include <cplib/fps/product_tree.hpp>
#include <cplib/fps/shift_of_sampling_points.hpp>
#include <cplib/fps/taylor_shift.hpp>
namespace cplib {
namespace detail {
// Bが2冪のとき(iB)!をサンプリング点シフトで一括構築。O(M(B+blockCount))。
template<Modint T> std::vector<T> factorialBlocks(Int blockSize,Int blockCount){std::vector<T> samples={T(1)};for(Int width=1;width<blockSize;width*=2){auto extended=samples;auto tail=shiftOfSamplingPoints(samples,T(width),3*width);extended.insert(extended.end(),tail.begin(),tail.end());std::vector<T> next(2*width);for(Int i=0;i<Int(next.size());++i)next[i]=extended[2*i]*extended[2*i+1]*T((2*i+1)*width);samples=std::move(next);}if(blockCount>blockSize){auto tail=shiftOfSamplingPoints(samples,T(blockSize),blockCount-blockSize);samples.insert(samples.end(),tail.begin(),tail.end());}std::vector<T> out(blockCount+1);out[0]=1;for(Int i=0;i<blockCount;++i)out[i+1]=out[i]*samples[i]*T((i+1)*blockSize);return out;}
}
template<Modint T> class LargeFactorial {
    std::uint32_t modulus;Int maxN,blockSize;std::vector<T> blocks;
public:
    // 素数法。前計算O(M(B+N/B))、保存領域O(1+N/B)。
    explicit LargeFactorial(Int limit=-1,Int width=1024):modulus(T::umod()),maxN(limit==-1?Int(modulus)-1:std::min(limit,Int(modulus)-1)),blockSize(width){assert(limit>=-1&&width>0&&std::has_single_bit(UInt(width))&&modulus>=2);while(blockSize>std::max(Int(1),maxN))blockSize/=2;blocks=detail::factorialBlocks<T>(blockSize,maxN/blockSize);}
    T fact(Int n)const{assert(!blocks.empty()&&modulus==T::umod()&&n>=0);if(n>=modulus)return T(0);assert(n<=maxN);Int index=n/blockSize;T out=blocks[index];for(Int i=index*blockSize+1;i<=n;++i)out*=T(i);return out;}
};
template<Modint T> auto initLargeFactorial(Int maxN=-1,Int blockSize=1024){return LargeFactorial<T>(maxN,blockSize);}
template<Modint T> T fact(const LargeFactorial<T>& table,Int n){return table.fact(n);}
// 問い合わせをソートし、余り区間は多点評価でまとめる元アルゴリズム。
// NTT使用時O(sqrt(p)log p+Q log Q+Q log³p)。
template<Modint T> std::vector<T> manyFactorials(std::span<const Int> ns){std::vector<T> out(ns.size());Int modulus=T::umod(),maxN=-1;for(Int n:ns){assert(n>=0);if(n<modulus)maxN=std::max(maxN,n);}if(maxN<0)return out;if(maxN<=1024){std::vector<T> fact(maxN+1);fact[0]=1;for(Int i=1;i<=maxN;++i)fact[i]=fact[i-1]*T(i);for(Int i=0;i<Int(ns.size());++i)if(ns[i]<modulus)out[i]=fact[ns[i]];return out;}Int blockSize=1;while(blockSize*blockSize<=maxN)blockSize*=2;auto blocks=detail::factorialBlocks<T>(blockSize,maxN/blockSize);std::vector<std::pair<Int,Int>> queries;queries.reserve(ns.size());for(Int i=0;i<Int(ns.size());++i)if(ns[i]<modulus){out[i]=blocks[ns[i]/blockSize];queries.emplace_back(ns[i],i);}std::sort(queries.begin(),queries.end());std::vector<T> polynomial={T(1),T(1)};for(Int width=1;width<blockSize;width*=2){if(width<=32){for(auto [n,index]:queries)if(n&width){Int start=n-n%(2*width);T product=1;for(Int i=1;i<=width;++i)product*=T(start+i);out[index]*=product;}}else{std::vector<T> points;Int previous=-1;for(auto [n,index]:queries)if(n&width){Int start=n-n%(2*width);if(start!=previous){points.push_back(T(start));previous=start;}}std::vector<T> values(points.size());for(Int first=0;first<Int(points.size());){Int last=std::min(first+width,Int(points.size()));auto evaluated=multipointEvaluation(polynomial,std::vector<T>(points.begin()+first,points.begin()+last));std::copy(evaluated.begin(),evaluated.end(),values.begin()+first);first=last;}Int pointIndex=-1;previous=-1;for(auto [n,index]:queries)if(n&width){Int start=n-n%(2*width);if(start!=previous){++pointIndex;previous=start;}out[index]*=values[pointIndex];}}if(2*width<blockSize)polynomial=polynomial*taylorShift(polynomial,T(width));}return out;}
}
