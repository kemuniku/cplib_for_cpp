#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/graph/topologicalsort.hpp>
#include <cplib/graph/maxflow.hpp>
namespace cplib {
// 二部グラフの単位容量最大流によるDAGの最小パス被覆。
template<class G> requires(DirectedGraph<G>&&UnWeightedGraph<G>) Int dag_minimum_path_cover(const G& g){assert(isDAG(g));auto flow=initMaxFlow(g.len*2+2);for(Int i=0;i<g.len;++i){for(auto [j,c]:g.to_and_cost(i))flow.add_edge(i,g.len+j,1);flow.add_edge(g.len*2,i,1);flow.add_edge(g.len+i,g.len*2+1,1);}return g.len-flow.flow(g.len*2,g.len*2+1);}
}
