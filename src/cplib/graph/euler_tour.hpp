#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
namespace detail {
// 次数を確認して非再帰Hierholzer法で辺番号を列挙する。O(V+E)。
template <GraphTypes G> std::vector<Int> euler_walk(const G &g, Int start, bool closed) {
    if constexpr (StaticGraphTypes<G>)
        g.static_graph_initialized_check();
    assert(start >= -1 && (start == -1 || start < g.len));
    if (!g.edge_count())
        return {};
    std::vector<Int> degree(g.len);
    for (auto e : g.edge_info) {
        ++degree[e.src];
        if constexpr (DirectedGraph<G>)
            --degree[e.dst];
        else
            ++degree[e.dst];
    }
    Int requiredStart = -1, endpoints = 0;
    for (Int v = 0; v < g.len; ++v) {
        if constexpr (DirectedGraph<G>) {
            if (degree[v] == 1) {
                if (requiredStart != -1)
                    return {};
                requiredStart = v;
            } else if (degree[v] == -1)
                ++endpoints;
            else if (degree[v])
                return {};
        } else if (degree[v] % 2) {
            if (requiredStart == -1)
                requiredStart = v;
            ++endpoints;
        }
    }
    if constexpr (DirectedGraph<G>) {
        if (endpoints != Int(requiredStart != -1))
            return {};
    } else if (endpoints != 0 && endpoints != 2)
        return {};
    if (closed && requiredStart != -1)
        return {};
    Int root = start;
    if (root == -1)
        root = requiredStart != -1 ? requiredStart : g.edge_info[0].src;
    else if (requiredStart != -1) {
        if constexpr (DirectedGraph<G>) {
            if (root != requiredStart)
                return {};
        } else if (degree[root] % 2 == 0)
            return {};
    }
    std::vector<Int> next(g.len), vertices{root}, incoming{-1}, result;
    std::vector<bool> used(g.edge_count());
    result.reserve(g.edge_count());
    while (!vertices.empty()) {
        Int v = vertices.back();
        auto adj = g.adjacency(v);
        if (next[v] == Int(adj.size())) {
            vertices.pop_back();
            Int id = incoming.back();
            incoming.pop_back();
            if (id != -1)
                result.push_back(id);
            continue;
        }
        auto e = adj[next[v]++];
        if (used[e.id])
            continue;
        used[e.id] = true;
        vertices.push_back(e.dst);
        incoming.push_back(e.id);
    }
    if (Int(result.size()) != g.edge_count())
        result.clear();
    else
        std::reverse(result.begin(), result.end());
    return result;
}
}

// 全辺を一度ずつ通り始点に戻る閉路を、通過順の辺番号列で返す。O(V+E) 時間・領域。
// start = -1 なら始点を自動選択する。指定した始点から存在しなければ空列。辺がない場合も空列。
// 辺番号は add_edge の戻り値。無向辺はどちら向きにも通れる。孤立点は無視する。
// 自己ループ・多重辺に対応し、重みは無視する。静的グラフは build 済みとする。グラフは変更しない。
template <GraphTypes G> auto euler_tour(const G &g, Int start = -1) {
    return detail::euler_walk(g, start, true);
}

// 全辺を一度ずつ通る経路を、通過順の辺番号列で返す。始点と終点は異なってもよい。O(V+E) 時間・領域。
// start = -1 なら始点を自動選択する。指定した始点から存在しなければ空列。辺がない場合も空列。
// 辺番号・対応グラフ・非破壊性は euler_tour と同じ。
template <GraphTypes G> auto euler_trail(const G &g, Int start = -1) {
    return detail::euler_walk(g, start, false);
}
}
