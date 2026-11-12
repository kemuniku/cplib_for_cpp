#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/iterator_views.hpp>
#include <stdexcept>
#include <ranges>

namespace cplib {
// 指定辺を除いた独立したグラフを返す。辺IDは再採番する。時間・領域 O(V+E)。
template <GraphTypes G> G remove_edge(const G &g, Int id) {
    if (id < 0 || id >= g.edge_count())
        throw std::invalid_argument("辺番号が範囲外です");
    G result(g.size());
    for (Int i = 0; i < g.edge_count(); ++i) {
        if (i == id)
            continue;
        const auto &e = g.edge_info[i];
        if constexpr (UnWeightedGraph<G>)
            result.add_edge(e.src, e.dst);
        else
            result.add_edge(e.src, e.dst, e.cost);
    }
    if constexpr (G::is_static)
        result.build();
    return result;
}
// 各辺を一本ずつ除いたグラフを遅延列挙する。ビュー取得 O(1)、一件の時間・領域 O(V+E)。
// 元グラフはビューの使用中に生存し、変更されないこと。
template <GraphTypes G> auto remove_each_edge(const G &g) {
    return std::views::iota(Int(0), g.edge_count()) |
           std::views::transform([&g](Int id) { return remove_edge(g, id); });
}
}
