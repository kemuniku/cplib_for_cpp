#pragma once
#include <cplib/common.hpp>
#include <cplib/collections/private/staticRMQ_kernel.hpp>
namespace cplib {
template<class T> class StaticRMQ {
    std::vector<std::vector<T>> table;
    std::vector<T> prefix_product,suffix_product,V;
    static constexpr bool simd=std::is_same_v<T,std::int32_t>||std::is_same_v<T,std::int64_t>;
public:
    StaticRMQ()=default;
    // 64要素ブロックの累積最小値とスパーステーブルを構築する。O(n+n/64 log n)。
    explicit StaticRMQ(std::span<const T> values):V(values.begin(),values.end()){
        Int n=V.size();if(!n)return;Int blocks=(n-1)/64+1,levels=std::bit_width(UInt(blocks));prefix_product.resize(n);suffix_product.resize(n);table.resize(levels);
        for(Int k=0;k<levels;++k)table[k].resize(blocks-(Int(1)<<k)+1);
        std::vector<T*> rows;for(auto& row:table)rows.push_back(row.data());
        if constexpr(simd)cplib::detail::static_rmq_kernel::build(V.data(),prefix_product.data(),suffix_product.data(),rows.data(),n);
        else cplib::detail::static_rmq_kernel::build_scalar(V.data(),prefix_product.data(),suffix_product.data(),rows.data(),n);
    }
    Int minLeft(Int r,const T& lower)const{assert(0<=r&&r<=Int(V.size()));if(!r)return 0;Int result=r-1;if(prefix_product[result]<lower){while(V[result]>=lower)--result;return result+1;}Int finish=result>>6;if(!finish)return 0;for(Int k=std::bit_width(UInt(finish));k-->0;){Int width=Int(1)<<k;if(width<=finish&&table[k][finish-width]>=lower)finish-=width;}if(!finish)return 0;result=(finish<<6)-1;while(V[result]>=lower)--result;return result+1;}
    Int maxRight(Int l,const T& lower)const{Int n=V.size();assert(0<=l&&l<=n);if(l==n)return l;Int result=l;if(suffix_product[l]<lower){while(V[result]>=lower)++result;return result;}Int first=(l>>6)+1,blocks=table[0].size();if(first==blocks)return n;for(Int k=std::bit_width(UInt(blocks-first));k-->0;){Int width=Int(1)<<k;if(first+width<=blocks&&table[k][first]>=lower)first+=width;}if(first==blocks)return n;result=first<<6;while(V[result]>=lower)++result;return result;}
    // 半開区間の最小値を返す。O(1)、短い区間は最大64要素を走査する。
    T query(Int l,Int r)const{
        assert(0<=l&&l<r&&r<=Int(V.size()));Int last=r-1,a=l/64,b=last/64;
        if(a==b){if constexpr(simd)if(r-l>=Int(32/sizeof(T)))return cplib::detail::static_rmq_kernel::scan_dispatch(V.data()+l,r-l);return *std::min_element(V.begin()+l,V.begin()+r);}
        T result=std::min(suffix_product[l],prefix_product[last]);
        if(a+1<b){Int k=std::bit_width(UInt(b-a-1))-1;result=std::min(result,std::min(table[k][a+1],table[k][b-(Int(1)<<k)]));}return result;
    }
};
template<class T> StaticRMQ<T> initRMQ(std::span<const T> values){return StaticRMQ<T>(values);}
template<class T> StaticRMQ<T> initRMQ(const std::vector<T>& values){return StaticRMQ<T>(values);}
template<class T> T query(const StaticRMQ<T>& rmq,Int l,Int r){return rmq.query(l,r);}
}

namespace cplib {template<class T> Int minLeft(const StaticRMQ<T>& s,Int r,T x){return s.minLeft(r,x);}template<class T> Int maxRight(const StaticRMQ<T>& s,Int l,T x){return s.maxRight(l,x);}}
