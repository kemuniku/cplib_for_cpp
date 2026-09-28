#pragma once
#include <cplib/graph/lowlink.hpp>
namespace cplib {
struct TwoEdgeConnectedComponents {std::vector<std::vector<Int>> groups;std::vector<Int> component;UnWeightedUnDirectedGraph forest;};
// lowlinkから橋で縮約森を構築する。O(V)。
inline TwoEdgeConnectedComponents initTwoEdgeConnectedComponents(const LowLink& ll){TwoEdgeConnectedComponents out;out.component.resize(ll.ord.size());for(Int v:ll.preorder){Int p=ll.parent[v];if(p==-1||ll.low[v]>ll.ord[p]){out.component[v]=out.groups.size();out.groups.push_back({v});}else{out.component[v]=out.component[p];out.groups[out.component[v]].push_back(v);}}out.forest=initUnWeightedUnDirectedGraph(out.groups.size());for(auto [u,v]:ll.bridges)out.forest.add_edge(out.component[u],out.component[v]);return out;}
template<UnDirectedGraph G> auto initTwoEdgeConnectedComponents(const G& g){return initTwoEdgeConnectedComponents(initLowLink(g));}
}
