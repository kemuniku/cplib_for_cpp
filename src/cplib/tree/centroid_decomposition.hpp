#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
struct CentroidDecomposition {
    Int root = -1;
    std::vector<Int> parent, depth;
    std::vector<std::vector<Int>> children;
};

struct CentroidDecompositionTree {
    UnWeightedDirectedGraph tree;
    Int root = -1;
    std::vector<Int> size, depth;
};

namespace detail {
// 反復探索で各成分の重心を選ぶ。O(n log n)時間・O(n)領域。
template <bool ReturnTree, UnDirectedGraph G>
auto centroid_decomposition_impl(const G &g, Int root) {
    Int n = g.len;
    using Result = std::conditional_t<ReturnTree, CentroidDecompositionTree, CentroidDecomposition>;
    Result result;
    if constexpr (ReturnTree) {
        result.tree = initUnWeightedDirectedGraph(n);
        result.size.resize(n);
    } else {
        result.parent.resize(n);
        result.children.resize(n);
    }
    result.depth.resize(n);
    if (!n)
        return result;
    assert(0 <= root && root < n);
    std::vector<bool> removed(n);
    std::vector<Int> parent(n), size(n), largest(n), seen(n), order;
    Int stamp = 0;
    std::vector<std::pair<Int, Int>> tasks{{root, -1}};
    while (!tasks.empty()) {
        auto [start, decompositionParent] = tasks.back();
        tasks.pop_back();
        ++stamp;
        order.clear();
        order.push_back(start);
        parent[start] = -1;
        seen[start] = stamp;
        for (std::size_t index = 0; index < order.size(); ++index) {
            Int u = order[index];
            size[u] = 1;
            largest[u] = 0;
            for (auto [v, c] : g.to_and_cost(u)) {
                if (removed[v] || v == parent[u])
                    continue;
                assert(seen[v] != stamp);
                seen[v] = stamp;
                parent[v] = u;
                order.push_back(v);
            }
        }
        if (decompositionParent == -1)
            assert(Int(order.size()) == n);
        Int centroid = -1;
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            Int u = *it;
            if (std::max(largest[u], Int(order.size()) - size[u]) <= Int(order.size()) / 2)
                centroid = u;
            Int p = parent[u];
            if (p != -1) {
                size[p] += size[u];
                largest[p] = std::max(largest[p], size[u]);
            }
        }
        if constexpr (ReturnTree)
            result.size[centroid] = order.size();
        else
            result.parent[centroid] = decompositionParent;
        if (decompositionParent == -1)
            result.root = centroid;
        else {
            result.depth[centroid] = result.depth[decompositionParent] + 1;
            if constexpr (ReturnTree)
                result.tree.add_edge(decompositionParent, centroid);
            else
                result.children[decompositionParent].push_back(centroid);
        }
        removed[centroid] = true;
        for (auto [v, c] : g.to_and_cost(centroid))
            if (!removed[v])
                tasks.emplace_back(v, centroid);
    }
    return result;
}
}

// 無向木の重心分解木を O(N log N) 時間・O(N) 領域で構築する。
// root は探索開始点で、分解木の根とは限らない。辺の重みは使わない。
// 空の木では root = -1、各配列は空。静的グラフは事前に build する。
// 返り値の頂点番号は元の木と共通。返り値のrootではparent[root] = -1、depth[root] = 0。
template <UnDirectedGraph G>
CentroidDecomposition initCentroidDecomposition(const G &g, Int root = 0) {
    return detail::centroid_decomposition_impl<false>(g, root);
}

// 重心分解木の有向グラフ・根・部分木サイズ・深さを O(N log N) 時間・O(N) 領域で構築する。
// 有向辺は分解木の親から子へ張り、頂点番号は元の木と共通。
// size[v] は重心分解木で v を根とする部分木の頂点数（v 自身を含む）。
// depth[v] は重心分解木の根からの深さで、根の深さは 0。
// root は探索開始点で、分解木の根とは限らない。辺の重みは使わない。
// 空の木では根は -1、size と depth は空。静的グラフは事前に build する。
template <UnDirectedGraph G>
CentroidDecompositionTree initCentroidDecompositionTree(const G &g, Int root = 0) {
    return detail::centroid_decomposition_impl<true>(g, root);
}
}
