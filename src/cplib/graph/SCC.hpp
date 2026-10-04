#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
// 強連結成分をトポロジカル順に返す。Kosaraju法、O(V+E)。
template <DirectedGraph G>
    requires UnWeightedGraph<G>
std::vector<std::vector<Int>> SCC(const G &g) {
    std::vector<Int> postorder(g.len, -1);
    std::vector<bool> used(g.len);
    Int count = g.len - 1;
    auto fdfs = [&](auto &&self, Int x) -> void {
        for (Int i : g[x])
            if (!used[i]) {
                used[i] = true;
                self(self, i);
            }
        postorder[count--] = x;
    };
    for (Int i = 0; i < g.len; ++i)
        if (!used[i]) {
            used[i] = true;
            fdfs(fdfs, i);
        }
    std::vector<std::vector<Int>> gout(g.len), groups;
    for (Int i = 0; i < g.len; ++i)
        for (Int j : g[i])
            gout[j].push_back(i);
    std::fill(used.begin(), used.end(), false);
    auto sdfs = [&](auto &&self, Int x) -> void {
        groups.back().push_back(x);
        for (Int i : gout[x])
            if (!used[i]) {
                used[i] = true;
                self(self, i);
            }
    };
    for (Int i : postorder)
        if (!used[i]) {
            used[i] = true;
            groups.emplace_back();
            sdfs(sdfs, i);
        }
    return groups;
}

// 縮約グラフ・頂点対応・成分一覧を返す。多重辺を保持する。O(V+E)。
// 強連結成分分解をします。
// 結果を、(頂点をまとめたグラフ,元の頂点→新頂点への対応,新頂点に含まれる頂点一覧)で返します。
template <DirectedGraph G>
    requires UnWeightedGraph<G>
auto SCCG(const G &g) {
    auto groups = SCC(g);
    std::vector<Int> mapping(g.len, -1);
    for (Int i = 0; i < Int(groups.size()); ++i)
        for (Int j : groups[i])
            mapping[j] = i;
    G result(groups.size());
    for (Int i = 0; i < g.len; ++i)
        for (Int j : g[i])
            if (mapping[i] != mapping[j])
                result.add_edge(mapping[i], mapping[j]);
    if constexpr (StaticGraphTypes<G>)
        result.build();
    return std::tuple{std::move(result), std::move(mapping), std::move(groups)};
}
}
