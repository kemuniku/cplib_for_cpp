#pragma once
#include <cplib/graph/maxflow.hpp>

namespace cplib {
template <class Cap> struct LowerBoundMaxFlowEdge {
    Int src, dst;
    Cap lower, upper, flow;
};

template <class Cap = Int> class LowerBoundMaxFlow {
    Int n;
    std::vector<LowerBoundMaxFlowEdge<Cap>> edges;
    bool solved = false;

public:
    // size頂点の下限流量つき最大流グラフを構築する。容量型の省略時はInt。容量型は符号付き整数とする。O(1)。
    explicit LowerBoundMaxFlow(Int size) : n(size) {
        static_assert(std::is_signed_v<Cap>);
        assert(n >= 0);
    }

    // 流量の下限lower、上限upperの有向辺を追加し、辺番号を返す。償却O(1)。
    Int add_edge(Int src, Int dst, Cap lower, Cap upper) {
        assert(0 <= src && src < n && 0 <= dst && dst < n && Cap(0) <= lower && lower <= upper);
        Int id = edges.size();
        edges.push_back({src, dst, lower, upper, Cap(0)});
        solved = false;
        return id;
    }

    // 下限を満たす非負のsrc-dst最大流量を返す。実現不可能なら-1。O(V^2(V+E))。
    // 呼び出すたびに最初から計算する。負の流量のみ実現可能な場合も-1を返す。
    // 各頂点の下限流量の収支の中間値、正の収支の総和、最大流量はCapに収まること。
    Cap flow(Int src, Int dst) {
        assert(0 <= src && src < n && 0 <= dst && dst < n && src != dst);
        solved = false;
        auto auxiliary = initMaxFlow<Cap>(n + 2);
        std::vector<Cap> balance(n);
        for (auto e : edges) {
            auxiliary.add_edge(e.src, e.dst, e.upper - e.lower);
            if (e.src != e.dst) {
                balance[e.src] -= e.lower;
                balance[e.dst] += e.lower;
            }
        }
        Int back = auxiliary.add_edge(dst, src, std::numeric_limits<Cap>::max());
        Cap required = 0;
        for (Int v = 0; v < n; ++v)
            if (balance[v] > Cap(0)) {
                auxiliary.add_edge(n, v, balance[v]);
                required += balance[v];
            } else if (balance[v] < Cap(0))
                auxiliary.add_edge(v, n + 1, -balance[v]);
        if (auxiliary.flow(n, n + 1, required) != required)
            return Cap(-1);
        Cap initial = auxiliary.get_edge(back).flow;
        auto residual = initMaxFlow<Cap>(n);
        for (Int i = 0; i < Int(edges.size()); ++i) {
            auto e = edges[i];
            Cap extra = auxiliary.get_edge(i).flow;
            residual.add_edge(e.src, e.dst, e.upper - e.lower - extra);
            residual.add_edge(e.dst, e.src, extra);
        }
        Cap result = initial + residual.flow(src, dst, std::numeric_limits<Cap>::max() - initial);
        for (Int i = 0; i < Int(edges.size()); ++i)
            edges[i].flow = edges[i].lower + auxiliary.get_edge(i).flow -
                            residual.get_edge(2 * i + 1).flow + residual.get_edge(2 * i).flow;
        solved = true;
        return result;
    }

    // flow成功後のi番目の辺の上下限と流量を返す。辺の追加後は再計算が必要。O(1)。
    auto get_edge(Int i) const {
        assert(solved);
        return edges[i];
    }

    // flow成功後の全辺の上下限と流量を追加順に返す。O(E)。
    auto get_edges() const {
        assert(solved);
        return edges;
    }
};

template <class Cap = Int> auto initLowerBoundMaxFlow(Int n, Cap = Cap(0)) {
    return LowerBoundMaxFlow<Cap>(n);
}

template <class Cap> Int add_edge(LowerBoundMaxFlow<Cap> &g, Int s, Int t, Cap l, Cap u) {
    return g.add_edge(s, t, l, u);
}

template <class Cap> auto get_edge(const LowerBoundMaxFlow<Cap> &g, Int i) {
    return g.get_edge(i);
}

template <class Cap> auto get_edges(const LowerBoundMaxFlow<Cap> &g) {
    return g.get_edges();
}

template <class Cap> Cap flow(LowerBoundMaxFlow<Cap> &g, Int s, Int t) {
    return g.flow(s, t);
}
}
