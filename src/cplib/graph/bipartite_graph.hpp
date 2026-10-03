#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
// DFSで二部グラフを判定する。O(V+E)。
template <UnDirectedGraph G> bool is_bipartite_graph(const G &g) {
    std::vector<Int> color(g.len, -1);
    auto dfs = [&](auto &&self, Int x) -> bool {
        for (auto [y, cost] : g.to_and_cost(x)) {
            if (color[y] == -1) {
                color[y] = color[x] ^ 1;
                if (!self(self, y))
                    return false;
            } else if (color[y] != (color[x] ^ 1))
                return false;
        }
        return true;
    };
    for (Int i = 0; i < g.len; ++i)
        if (color[i] == -1) {
            color[i] = 0;
            if (!dfs(dfs, i))
                return false;
        }
    return true;
}
}
