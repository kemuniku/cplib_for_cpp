#pragma once
#include <cplib/graph/biconnected_components.hpp>

namespace cplib {
struct BlockCutTree {
    UnWeightedUnDirectedGraph forest;
    std::vector<Int> id, articulation;
    std::vector<std::vector<Int>> groups;
};

// 成分ノード[0, groups.size())と関節点ノードからなるblock-cut forestをO(V)で構築する。
// id[v]は関節点なら専用ノード、それ以外なら所属成分を指す。
// 関節点ノードgroups.size()+iはarticulation[i]に対応する。
inline BlockCutTree initBlockCutTree(const BiconnectedComponents &bc) {
    BlockCutTree out;
    out.groups = bc.groups;
    out.articulation = bc.articulation;
    out.id.resize(bc.belong.size());
    out.forest = initUnWeightedUnDirectedGraph(bc.groups.size() + bc.articulation.size());
    for (Int v = 0; v < Int(bc.belong.size()); ++v)
        out.id[v] = bc.belong[v][0];
    for (Int i = 0; i < Int(bc.articulation.size()); ++i) {
        Int v = bc.articulation[i], node = bc.groups.size() + i;
        out.id[v] = node;
        for (Int group : bc.belong[v])
            out.forest.add_edge(node, group);
    }
    return out;
}

// 無向グラフのblock-cut forestをO(V+E)時間・領域で構築します。
template <UnDirectedGraph G> auto initBlockCutTree(const G &g) {
    return initBlockCutTree(initBiconnectedComponents(g));
}
}
