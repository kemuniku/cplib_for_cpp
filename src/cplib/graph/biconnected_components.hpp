#pragma once
#include <cplib/graph/lowlink.hpp>

namespace cplib {
struct BiconnectedComponents {
    std::vector<std::vector<Int>> groups, belong;
    std::vector<Int> articulation;

    // 頂点vと接続する辺を削除した後の連結成分数から削除前の個数を引いた差分をO(1)で返します。
    Int component_count_delta_after_removal(Int v) const {
        return groups[belong[v][0]].size() == 1 ? -1 : Int(belong[v].size()) - 1;
    }
};

inline Int component_count_delta_after_removal(const BiconnectedComponents &bc, Int v) {
    return bc.component_count_delta_after_removal(v);
}

// 計算済みlowlinkの帰りがけ順から成分をまとめる。O(V)。
// 計算済みlowlinkから二重頂点連結成分の頂点集合をO(V)で求めます。
// 橋の両端も一成分とし、自己ループは無視します。孤立点は単独の成分です。
inline BiconnectedComponents initBiconnectedComponents(const LowLink &ll) {
    BiconnectedComponents out;
    out.articulation = ll.articulation;
    out.belong.resize(ll.ord.size());
    std::vector<Int> pending;
    for (Int v : ll.postorder) {
        Int p = ll.parent[v];
        std::vector<Int> group;
        if (p == -1) {
            if (out.belong[v].empty())
                group.push_back(v);
        } else if (ll.low[v] >= ll.ord[p]) {
            while (!pending.empty() && ll.ord[pending.back()] > ll.ord[v]) {
                group.push_back(pending.back());
                pending.pop_back();
            }
            group.push_back(v);
            group.push_back(p);
        } else
            pending.push_back(v);
        if (!group.empty()) {
            Int id = out.groups.size();
            for (Int u : group)
                out.belong[u].push_back(id);
            out.groups.push_back(std::move(group));
        }
    }
    return out;
}

// グラフから直接分解する非再帰版。親辺もlow更新に含める元実装を保持。O(V+E)。
// 無向グラフをO(V+E)時間・領域で二重頂点連結成分に分解します。
// 静的グラフはbuild済みとします。自己ループは無視し、DFSは非再帰です。
template <UnDirectedGraph G> BiconnectedComponents initBiconnectedComponents(const G &g) {
    if constexpr (StaticGraphTypes<G>)
        g.static_graph_initialized_check();
    Int n = g.len, timer = 0;
    BiconnectedComponents out;
    out.belong.resize(n);
    std::vector<Int> ord(n, -1), low(n), parent(n, -1), next(n), pending;
    for (Int root = 0; root < n; ++root) {
        if (ord[root] != -1)
            continue;
        Int v = root;
        ord[v] = low[v] = timer++;
        while (v != -1) {
            auto adj = g.adjacency(v);
            if (next[v] < Int(adj.size())) {
                Int to = adj[next[v]++].dst;
                if (ord[to] == -1) {
                    parent[to] = v;
                    ord[to] = low[to] = timer++;
                    v = to;
                } else
                    low[v] = std::min(low[v], ord[to]);
            } else {
                Int p = parent[v];
                std::vector<Int> group;
                if (p == -1) {
                    if (out.belong[v].empty())
                        group.push_back(v);
                } else {
                    low[p] = std::min(low[p], low[v]);
                    if (low[v] >= ord[p]) {
                        Int first = pending.size();
                        while (first > 0 && ord[pending[first - 1]] > ord[v])
                            --first;
                        Int size = pending.size() - first;
                        group.resize(size + 2);
                        for (Int i = 0; i < size; ++i)
                            group[i] = pending[pending.size() - 1 - i];
                        pending.resize(first);
                        group[size] = v;
                        group[size + 1] = p;
                    } else
                        pending.push_back(v);
                }
                if (!group.empty()) {
                    Int id = out.groups.size();
                    for (Int u : group)
                        out.belong[u].push_back(id);
                    out.groups.push_back(std::move(group));
                }
                v = p;
            }
        }
    }
    for (Int v = 0; v < n; ++v)
        if (out.belong[v].size() > 1)
            out.articulation.push_back(v);
    return out;
}
}
