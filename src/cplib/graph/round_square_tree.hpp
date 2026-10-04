#pragma once
#include <cplib/graph/biconnected_components.hpp>

namespace cplib {
// 元頂点[0,n)と成分ノード[n,n+groups.size())を結ぶ円方木（森）。O(V)。
inline UnWeightedUnDirectedGraph initRoundSquareTree(const BiconnectedComponents &bc) {
    Int n = bc.belong.size();
    auto out = initUnWeightedUnDirectedGraph(n + bc.groups.size());
    for (Int i = 0; i < Int(bc.groups.size()); ++i)
        for (Int v : bc.groups[i])
            out.add_edge(v, n + i);
    return out;
}

// 無向グラフの円方木をO(V+E)時間・領域で構築する。自己ループと重みは無視する。
// 静的グラフはbuild済みとする。成分分解は非再帰。
template <UnDirectedGraph G> auto initRoundSquareTree(const G &g) {
    return initRoundSquareTree(initBiconnectedComponents(g));
}
}
