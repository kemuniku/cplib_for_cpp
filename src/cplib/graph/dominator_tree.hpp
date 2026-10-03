#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
// 各頂点の直近支配頂点を返す。根は自身、根から到達不能な頂点は-1。根は有効な頂点とする。
// Lengauer-Tarjan法でO((V+E)log V)時間、O(V+E)領域。DFS・経路圧縮ともに非再帰。
// 自己ループ・多重辺に対応する。重みは無視し、静的グラフはbuild済みとする。
template <class G> std::vector<Int> dominator_tree(const G &g, Int root) {
    Int n;
    if constexpr (GraphTypes<G>)
        n = g.len;
    else
        n = g.size();
    assert(0 <= root && root < n);
    auto degree = [&](Int u) -> Int {
        if constexpr (GraphTypes<G>)
            return g.adjacency(u).size();
        else
            return g[u].size();
    };
    auto neighbor = [&](Int u, Int i) -> Int {
        if constexpr (GraphTypes<G>)
            return g.adjacency(u)[i].dst;
        else
            return g[u][i];
    };
    std::vector<Int> order(n, -1), vertex, parent(n), next(n), stack;
    std::vector<std::vector<Int>> pred(n);
    order[root] = 0;
    vertex.push_back(root);
    stack.push_back(root);
    while (!stack.empty()) {
        Int u = stack.back();
        if (next[u] == degree(u)) {
            stack.pop_back();
            continue;
        }
        Int v = neighbor(u, next[u]++);
        assert(0 <= v && v < n);
        pred[v].push_back(u);
        if (order[v] < 0) {
            order[v] = vertex.size();
            parent[vertex.size()] = order[u];
            vertex.push_back(v);
            stack.push_back(v);
        }
    }
    Int m = vertex.size();
    std::vector<Int> semi(m), label(m), ancestor(m, -1), idom(m), head(m, -1), link(m);
    std::iota(semi.begin(), semi.end(), 0);
    label = semi;
    auto eval = [&](Int v) {
        Int x = v;
        while (ancestor[x] >= 0 && ancestor[ancestor[x]] >= 0) {
            stack.push_back(x);
            x = ancestor[x];
        }
        while (!stack.empty()) {
            Int u = stack.back();
            stack.pop_back();
            Int p = ancestor[u];
            if (semi[label[p]] < semi[label[u]])
                label[u] = label[p];
            ancestor[u] = ancestor[p];
        }
        return label[v];
    };
    for (Int w = m - 1; w > 0; --w) {
        for (Int v : pred[vertex[w]])
            semi[w] = std::min(semi[w], semi[eval(order[v])]);
        link[w] = head[semi[w]];
        head[semi[w]] = w;
        Int p = parent[w];
        ancestor[w] = p;
        for (Int v = head[p]; v != -1; v = link[v]) {
            Int u = eval(v);
            idom[v] = semi[u] < semi[v] ? u : p;
        }
        head[p] = -1;
    }
    for (Int w = 1; w < m; ++w)
        if (idom[w] != semi[w])
            idom[w] = idom[idom[w]];
    std::vector<Int> out(n, -1);
    out[root] = root;
    for (Int w = 1; w < m; ++w)
        out[vertex[w]] = vertex[idom[w]];
    return out;
}
}
