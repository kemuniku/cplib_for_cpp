#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
// 全方位木DP。元実装の結合順を保持。演算O(1)なら時間・空間O(N)。
// Eはモノイドの型、put_vertexの戻り値は頂点DPの型V。
// mergeはモノイドの結合、eは単位元、put_edgeは辺(u, v)の情報、put_vertexは頂点vの情報を付与する。
template <class E, class G, class Merge, class PutEdge, class PutVertex>
    requires(UnDirectedGraph<G> && UnWeightedGraph<G>)
std::vector<E> solve_Rerooting_raw(const G &g, Merge merge, E e, PutEdge put_edge,
                                   PutVertex put_vertex) {
    std::vector<std::vector<E>> L(g.len), R(g.len);
    std::vector<E> result(g.len);
    if (g.len == 0)
        return result;
    auto dfs1 = [&](auto &&self, Int x, Int p) -> E {
        std::vector<E> values{e};
        for (auto [y, c] : g.to_and_cost(x))
            if (y != p)
                values.push_back(put_edge(put_vertex(self(self, y, x), y), x, y));
        values.push_back(e);
        E now = e;
        auto &l = L[x];
        auto &r = R[x];
        l.resize(values.size());
        r.resize(values.size());
        for (Int i = 0; i < Int(values.size()); ++i) {
            now = merge(now, values[i]);
            l[i] = now;
        }
        now = e;
        for (Int i = Int(values.size()) - 1; i >= 0; --i) {
            now = merge(values[i], now);
            r[i] = now;
        }
        return now;
    };
    result[0] = dfs1(dfs1, 0, -1);
    auto dfs2 = [&](auto &&self, Int x, Int p, E value) -> void {
        if (p != -1)
            result[x] = merge(L[x].back(), value);
        Int i = 0;
        for (auto [y, c] : g.to_and_cost(x))
            if (y != p) {
                E a = merge(L[x][i], R[x][i + 2]);
                if (p != -1)
                    a = merge(a, value);
                self(self, y, x, put_edge(put_vertex(a, x), y, x));
                ++i;
            }
    };
    dfs2(dfs2, 0, -1, e);
    return result;
}

// 各頂点を根とするDP値Vを返す。
// mergeはモノイドEの結合、eは単位元、put_edgeは辺の情報、put_vertexは頂点の情報を付与する。
template <class E, class G, class Merge, class PutEdge, class PutVertex>
    requires(UnDirectedGraph<G> && UnWeightedGraph<G>)
auto solve_Rerooting(const G &g, Merge merge, E e, PutEdge put_edge, PutVertex put_vertex) {
    auto raw = solve_Rerooting_raw(g, merge, e, put_edge, put_vertex);
    using V = std::decay_t<std::invoke_result_t<PutVertex, E, Int>>;
    std::vector<V> result;
    result.reserve(raw.size());
    for (Int i = 0; i < Int(raw.size()); ++i)
        result.push_back(put_vertex(raw[i], i));
    return result;
}
}
