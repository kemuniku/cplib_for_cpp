#pragma once
#include <cplib/graph/graph.hpp>
namespace cplib {
// 辺番号とグラフ形式を維持して逆向きにする。O(V+E)。
template<DirectedGraph G> G reverse_edge(const G& g){G result(g.len);for(auto e:g.edge_info){if constexpr(UnWeightedGraph<G>)result.add_edge(e.dst,e.src);else result.add_edge(e.dst,e.src,e.cost);}if constexpr(StaticGraphTypes<G>)result.build();return result;}
}
