#pragma once
#include <cplib/graph/graph.hpp>
#include <queue>
namespace cplib {
// Prüfer列から木を復元する。次数と頂点番号のヒープを用いO(n log n)。
inline UnWeightedUnDirectedGraph prufer_decode(std::span<const Int> a){
    Int n=a.size()+2;auto result=initUnWeightedUnDirectedGraph(n);std::vector<Int> cnt(n,1);for(Int v:a){assert(0<=v&&v<n);++cnt[v];}
    std::priority_queue<std::pair<Int,Int>,std::vector<std::pair<Int,Int>>,std::greater<>> q;for(Int i=0;i<n;++i)q.emplace(cnt[i],i);
    for(Int v:a){auto [c,u]=q.top();q.pop();while(cnt[u]!=c){std::tie(c,u)=q.top();q.pop();}result.add_edge(u,v);--cnt[u];--cnt[v];if(cnt[u])q.emplace(cnt[u],u);if(cnt[v])q.emplace(cnt[v],v);}
    std::vector<Int> u;for(Int i=0;i<n;++i)if(cnt[i]==1)u.push_back(i);result.add_edge(u[0],u[1]);return result;
}
}
