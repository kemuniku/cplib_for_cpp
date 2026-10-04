#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
// 通過順の閉路の辺番号列から頂点列を復元する。O(L) 時間・領域。空列には空列を返す。
// 入力は cycle_detection の結果と同様に頂点・辺が重複しない閉路とする。
// edges[i] は結果の i 番目から (i+1) mod L 番目への辺。末尾に始点を重複させない。
template <GraphTypes G>
std::vector<Int> restore_cycle_vertices(const G &g, std::span<const Int> edges) {
    if (edges.empty())
        return {};
    auto first = g.edge_info[edges[0]];
    Int v = first.src;
    if constexpr (UnDirectedGraph<G>) {
        if (edges.size() > 1) {
            auto second = g.edge_info[edges[1]];
            if (first.dst != second.src && first.dst != second.dst)
                v = first.dst;
        }
    }
    std::vector<Int> out;
    out.reserve(edges.size());
    for (Int id : edges) {
        out.push_back(v);
        auto e = g.edge_info[id];
        if constexpr (DirectedGraph<G>)
            v = e.dst;
        else
            v = e.src == v ? e.dst : e.src;
    }
    return out;
}

// 閉路を一つ、通過順の辺番号列で返す。存在しなければ空列。O(V+E) 時間、O(V) 追加領域。
// 閉路内の頂点・辺は重複しない。各辺の端点は get_edge で取得できる。
// 自己ループ・多重辺に対応する非再帰 DFS。重みは無視し、静的グラフは build 済みとする。
template <GraphTypes G> std::vector<Int> cycle_detection(const G &g) {
    if constexpr (StaticGraphTypes<G>)
        g.static_graph_initialized_check();
    std::vector<std::uint8_t> state(g.len);
    std::vector<Int> next(g.len), parentEdge(g.len, -1), position(g.len), stack;
    for (Int root = 0; root < g.len; ++root) {
        if (state[root])
            continue;
        state[root] = 1;
        stack.push_back(root);
        position[root] = 0;
        while (!stack.empty()) {
            Int v = stack.back();
            auto adj = g.adjacency(v);
            if (next[v] == Int(adj.size())) {
                state[v] = 2;
                stack.pop_back();
                continue;
            }
            auto e = adj[next[v]++];
            Int to = e.dst, id = e.id;
            if constexpr (UnDirectedGraph<G>)
                if (id == parentEdge[v])
                    continue;
            if (state[to] == 0) {
                parentEdge[to] = id;
                position[to] = stack.size();
                state[to] = 1;
                stack.push_back(to);
            } else if (state[to] == 1) {
                std::vector<Int> out;
                for (Int i = position[to] + 1; i < Int(stack.size()); ++i)
                    out.push_back(parentEdge[stack[i]]);
                out.push_back(id);
                return out;
            }
        }
    }
    return {};
}
}
