#pragma once
#include <cplib/graph/graph.hpp>
namespace cplib {
// 木を二回DFSして直径長と端点を返す。O(V)。元実装と同じく非負辺を前提とする。
template<UnDirectedGraph G> auto diameter_and_edge(const G& g){
    using Cost=typename G::cost_type;Int u=0,v=0;Cost cur=0;
    auto dfs=[&](auto&& self,Int x,Int par,Cost d,Int& end)->void{if(d>cur){cur=d;end=x;}for(auto [y,cost]:g.to_and_cost(x))if(y!=par)self(self,y,x,d+cost,end);};
    dfs(dfs,0,-1,Cost(0),u);cur=0;dfs(dfs,u,-1,Cost(0),v);return std::tuple{cur,u,v};
}
template<UnDirectedGraph G> auto diameter(const G& g){return std::get<0>(diameter_and_edge(g));}
// 直径上の頂点列を復元する。O(V)。
template<UnDirectedGraph G> auto diameter_path(const G& g){auto [d,u,v]=diameter_and_edge(g);std::vector<Int> path;auto dfs=[&,target=v](auto&& self,Int x,Int par)->bool{path.push_back(x);if(x==target)return true;for(auto [y,cost]:g.to_and_cost(x))if(y!=par&&self(self,y,x))return true;path.pop_back();return false;};dfs(dfs,u,-1);return std::pair{d,path};}
}
