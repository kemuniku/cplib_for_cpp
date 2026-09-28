#pragma once
#include <cplib/graph/private/warshall_floyd_common.hpp>
#include <cplib/collections/unionfind.hpp>
namespace cplib {
namespace detail {
template<WeightedGraph G> auto kruskal_edges(const G& g){std::vector<std::tuple<typename G::cost_type,Int,Int>> edges;for(Int i=0;i<g.len;++i)for(auto [j,c]:g[i])if(i<j)edges.emplace_back(c,i,j);std::sort(edges.begin(),edges.end());return edges;}
}
// 最小全域木の重みを求める。非連結ならinf。O(E log E)。
template<WeightedGraph G,class T=typename G::cost_type> requires UnDirectedGraph<G> T get_MST_cost(const G& g,T zero=T(0),T inf=detail::distance_inf<T>()){
    auto edges=detail::kruskal_edges(g);UnionFind uf(g.len);T result=zero;for(auto [c,i,j]:edges)if(!uf.issame(i,j)){result+=c;uf.unite(i,j);}return uf.count==1?result:inf;
}
// 最小全域森を動的無向グラフで返す。O(E log E)。
template<WeightedGraph G> requires UnDirectedGraph<G> auto to_MST_Graph(const G& g){auto result=initWeightedUnDirectedGraph<typename G::cost_type>(g.len);UnionFind uf(g.len);for(auto [c,i,j]:detail::kruskal_edges(g))if(!uf.issame(i,j)){result.add_edge(i,j,c);uf.unite(i,j);}return result;}
}
