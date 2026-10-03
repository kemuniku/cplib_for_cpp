#pragma once
#include <cplib/graph/private/flow_graph.hpp>

namespace cplib {
template <class Cap> using MaxFlowEdge = FlowEdge<Cap>;

template <class Cap = Int> class MaxFlow : public detail::FlowGraph<Cap> {
    using Base = detail::FlowGraph<Cap>;
    using Base::graph;

public:
    explicit MaxFlow(Int n) : Base(n) {
    }

    // Dinic。終点からの非再帰DFSで流量をまとめて返す。O(V²E)。
    // Dinic法でlimit以下の流量を追加し、追加流量を返す。O(V^2 E)。
    Cap flow(Int src, Int dst, Cap limit = std::numeric_limits<Cap>::max()) {
        Int n = graph.size();
        assert(0 <= src && src < n && 0 <= dst && dst < n && src != dst && limit >= Cap(0));
        Cap result = 0;
        std::vector<Int> level(n), iter(n), queue(n);
        std::vector<Cap> requested(n), sent(n);
        while (result < limit) {
            std::fill(level.begin(), level.end(), -1);
            std::fill(iter.begin(), iter.end(), 0);
            level[src] = 0;
            queue[0] = src;
            Int head = 0, tail = 1;
            bool reached = false;
            while (head < tail && !reached) {
                Int v = queue[head++];
                for (auto e : graph[v])
                    if (e.cap > Cap(0) && level[e.dst] < 0) {
                        level[e.dst] = level[v] + 1;
                        if (e.dst == dst) {
                            reached = true;
                            break;
                        }
                        queue[tail++] = e.dst;
                    }
            }
            if (level[dst] < 0)
                break;
            Int depth = 0;
            queue[0] = dst;
            requested[0] = limit - result;
            sent[0] = 0;
            while (depth >= 0) {
                Int v = queue[depth];
                if (v == src)
                    sent[depth] = requested[depth];
                else {
                    while (iter[v] < Int(graph[v].size()) && sent[depth] < requested[depth]) {
                        auto e = graph[v][iter[v]];
                        if (level[e.dst] >= 0 && level[e.dst] < level[v] &&
                            graph[e.dst][e.rev].cap > Cap(0))
                            break;
                        ++iter[v];
                    }
                    if (sent[depth] < requested[depth] && iter[v] < Int(graph[v].size())) {
                        auto e = graph[v][iter[v]];
                        requested[depth + 1] =
                            std::min(requested[depth] - sent[depth], graph[e.dst][e.rev].cap);
                        sent[depth + 1] = 0;
                        queue[++depth] = e.dst;
                        continue;
                    }
                    if (sent[depth] < requested[depth])
                        level[v] = n;
                }
                Cap pushed = sent[depth--];
                if (depth < 0) {
                    result += pushed;
                    break;
                }
                Int parent = queue[depth], i = iter[parent], rev = graph[parent][i].rev;
                graph[parent][i].cap += pushed;
                graph[v][rev].cap -= pushed;
                sent[depth] += pushed;
                if (sent[depth] < requested[depth])
                    ++iter[parent];
            }
        }
        return result;
    }
};

template <class Cap = Int> auto initMaxFlow(Int n, Cap = Cap(0)) {
    return MaxFlow<Cap>(n);
}

template <class Cap> Int add_edge(MaxFlow<Cap> &g, Int s, Int t, Cap c) {
    return g.add_edge(s, t, c);
}

template <class Cap> auto get_edge(const MaxFlow<Cap> &g, Int i) {
    return g.get_edge(i);
}

template <class Cap> auto get_edges(const MaxFlow<Cap> &g) {
    return g.get_edges();
}

template <class Cap>
Cap flow(MaxFlow<Cap> &g, Int s, Int t, Cap limit = std::numeric_limits<Cap>::max()) {
    return g.flow(s, t, limit);
}

// 残余グラフでsから到達可能な頂点を返す。最大流計算後は最小カット。O(V+E)。
template <class Cap> auto min_cut(const MaxFlow<Cap> &g, Int s) {
    return g.min_cut(s);
}
}
