#pragma once
#include <cplib/graph/graph.hpp>
namespace cplib {
namespace detail {
// 次数を確認して非再帰Hierholzer法で辺番号を列挙する。O(V+E)。
template<GraphTypes G> std::vector<Int> euler_walk(const G& g,Int start,bool closed){
    if constexpr(StaticGraphTypes<G>)g.static_graph_initialized_check();
    assert(start>=-1 && (start==-1 || start<g.len));if(!g.edge_count())return {};
    std::vector<Int> degree(g.len);
    for(auto e:g.edge_info){++degree[e.src];if constexpr(DirectedGraph<G>)--degree[e.dst];else ++degree[e.dst];}
    Int requiredStart=-1,endpoints=0;
    for(Int v=0;v<g.len;++v){
        if constexpr(DirectedGraph<G>){if(degree[v]==1){if(requiredStart!=-1)return {};requiredStart=v;}else if(degree[v]==-1)++endpoints;else if(degree[v])return {};}
        else if(degree[v]%2){if(requiredStart==-1)requiredStart=v;++endpoints;}
    }
    if constexpr(DirectedGraph<G>){if(endpoints!=Int(requiredStart!=-1))return {};}
    else if(endpoints!=0 && endpoints!=2)return {};
    if(closed && requiredStart!=-1)return {};
    Int root=start;
    if(root==-1)root=requiredStart!=-1?requiredStart:g.edge_info[0].src;
    else if(requiredStart!=-1){if constexpr(DirectedGraph<G>){if(root!=requiredStart)return {};}else if(degree[root]%2==0)return {};}
    std::vector<Int> next(g.len),vertices{root},incoming{-1},result;std::vector<bool> used(g.edge_count());result.reserve(g.edge_count());
    while(!vertices.empty()){
        Int v=vertices.back();auto adj=g.adjacency(v);
        if(next[v]==Int(adj.size())){vertices.pop_back();Int id=incoming.back();incoming.pop_back();if(id!=-1)result.push_back(id);continue;}
        auto e=adj[next[v]++];if(used[e.id])continue;used[e.id]=true;vertices.push_back(e.dst);incoming.push_back(e.id);
    }
    if(Int(result.size())!=g.edge_count())result.clear();else std::reverse(result.begin(),result.end());return result;
}
}
// 全辺を一度ずつ通る閉路の辺番号列を返す。存在しなければ空列。O(V+E)。
template<GraphTypes G> auto euler_tour(const G& g,Int start=-1){return detail::euler_walk(g,start,true);}
// 全辺を一度ずつ通る経路の辺番号列を返す。存在しなければ空列。O(V+E)。
template<GraphTypes G> auto euler_trail(const G& g,Int start=-1){return detail::euler_walk(g,start,false);}
}
