#pragma once
#include <cplib/graph/graph.hpp>
#include <iostream>
#include <sstream>

namespace cplib {
namespace detail {
template <GraphTypes G, class F> void graph_debug_edges(const G &g, F f) {
    for (Int x = 0; x < g.len; ++x)
        for (auto [y, c] : g.to_and_cost(x))
            if (G::is_directed || y >= x)
                f(x, y, c);
}
}

// 元実装と同じ隣接順で出力。無向自己ループは隣接配列どおり2回出力。O(V+E)。
// 頂点番号に indexed を加えてグラフを出力する。O(V + E)。
template <GraphTypes G>
void dump_graph(const G &g, Int indexed = 0, std::ostream &output = std::cout) {
    Int m = 0;
    detail::graph_debug_edges(g, [&](Int, Int, auto) { ++m; });
    output << g.len << ' ' << m << '\n';
    detail::graph_debug_edges(g, [&](Int x, Int y, auto c) {
        output << x + indexed << ' ' << y + indexed;
        if constexpr (WeightedGraph<G>)
            output << ' ' << c;
        output << '\n';
    });
}

template <GraphTypes G> void dump_graph(const G &g, std::ostream &output) {
    dump_graph(g, 0, output);
}

template <GraphTypes G> std::string to_graph_graph(const G &g, bool indexed = false) {
    Int m = 0;
    detail::graph_debug_edges(g, [&](Int, Int, auto) { ++m; });
    std::ostringstream out;
    out << "https://hello-world-494ec.firebaseapp.com/?format=normal&indexed="
        << (indexed ? "true" : "false") << "&weighted=" << (WeightedGraph<G> ? "true" : "false")
        << "&directed=" << (DirectedGraph<G> ? "true" : "false") << "&data=" << g.len << '+' << m;
    detail::graph_debug_edges(g, [&](Int x, Int y, auto c) {
        out << "%0A" << x + Int(indexed) << '+' << y + Int(indexed);
        if constexpr (WeightedGraph<G>)
            out << '+' << c;
    });
    return out.str();
}
}
