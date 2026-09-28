#pragma once
#include <cplib/collections/deletable_heapqueue.hpp>
namespace cplib {
class TopK_sum_heapq {
    Deletable_HeapQueue<Int> G,L;
public:
    Int sm=0,topk=0,k;
    TopK_sum_heapq(std::span<const Int> input,Int k):k(k){assert(0<=k&&k<=Int(input.size()));std::vector<Int> v(input.begin(),input.end());std::sort(v.begin(),v.end(),std::greater<Int>());for(Int i=0;i<k;++i){G.push(v[i]);sm+=v[i];topk+=v[i];}for(Int i=k;i<Int(v.size());++i){L.push(-v[i]);sm+=v[i];}}
    // 大きいk個と残りを2つの削除可能ヒープで保持。償却O(log N)。
    void incl(Int x){sm+=x;if(G.len()<k){topk+=x;G.push(x);}else if(k==0)L.push(-x);else if(G[0]<=x){topk+=x;G.push(x);Int tmp=G.pop();L.push(-tmp);topk-=tmp;}else L.push(-x);}
    void excl(Int x){sm-=x;if(k==0)L.erase(-x);else if(G[0]<=x){topk-=x;G.erase(x);if(L.len()>0){Int tmp=-L.pop();G.push(tmp);topk+=tmp;}}else L.erase(-x);}
    void addK(){Int tmp=-L.pop();topk+=tmp;++k;G.push(tmp);}void minusK(){Int tmp=G.pop();topk-=tmp;--k;L.push(-tmp);}
    // 変更幅に比例する回数のヒープ操作。
    void setK(Int value){while(k>value)minusK();while(k<value)addK();}
};
inline auto initTopKHeapq(std::span<const Int> v,Int k){return TopK_sum_heapq(v,k);}
inline void incl(TopK_sum_heapq& q,Int x){q.incl(x);}inline void excl(TopK_sum_heapq& q,Int x){q.excl(x);}inline void addK(TopK_sum_heapq& q){q.addK();}inline void minusK(TopK_sum_heapq& q){q.minusK();}inline void setK(TopK_sum_heapq& q,Int k){q.setK(k);}
}
