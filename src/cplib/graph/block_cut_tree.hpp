#pragma once
#include <cplib/graph/biconnected_components.hpp>
namespace cplib {
struct BlockCutTree {UnWeightedUnDirectedGraph forest;std::vector<Int> id,articulation;std::vector<std::vector<Int>> groups;};
// 成分ノードの後ろに関節点専用ノードを追加する。O(V)。
inline BlockCutTree initBlockCutTree(const BiconnectedComponents& bc){BlockCutTree out;out.groups=bc.groups;out.articulation=bc.articulation;out.id.resize(bc.belong.size());out.forest=initUnWeightedUnDirectedGraph(bc.groups.size()+bc.articulation.size());for(Int v=0;v<Int(bc.belong.size());++v)out.id[v]=bc.belong[v][0];for(Int i=0;i<Int(bc.articulation.size());++i){Int v=bc.articulation[i],node=bc.groups.size()+i;out.id[v]=node;for(Int group:bc.belong[v])out.forest.add_edge(node,group);}return out;}
template<UnDirectedGraph G> auto initBlockCutTree(const G& g){return initBlockCutTree(initBiconnectedComponents(g));}
}
