#pragma once
#include <cplib/common.hpp>

namespace cplib {
template <class Cap> struct FlowEdge {
    Int src, dst;
    Cap cap, flow;
};

namespace detail {
template <class Cap> class FlowGraph {
protected:
    struct Arc {
        Int dst, rev;
        Cap cap;
    };

    std::vector<std::vector<Arc>> graph;
    std::vector<std::pair<Int, Int>> positions;

public:
    explicit FlowGraph(Int n) : graph(n) {
        assert(n >= 0);
    }

    Int add_edge(Int src, Int dst, Cap cap) {
        assert(0 <= src && src < Int(graph.size()) && 0 <= dst && dst < Int(graph.size()) &&
               cap >= Cap(0));
        Int id = positions.size(), index = graph[src].size(),
            rev = graph[dst].size() + Int(src == dst);
        positions.emplace_back(src, index);
        graph[src].push_back({dst, rev, cap});
        graph[dst].push_back({src, index, Cap(0)});
        return id;
    }

    FlowEdge<Cap> get_edge(Int i) const {
        auto [src, index] = positions[i];
        auto e = graph[src][index];
        Cap flow = graph[e.dst][e.rev].cap;
        return {src, e.dst, e.cap + flow, flow};
    }

    std::vector<FlowEdge<Cap>> get_edges() const {
        std::vector<FlowEdge<Cap>> result;
        result.reserve(positions.size());
        for (Int i = 0; i < Int(positions.size()); ++i)
            result.push_back(get_edge(i));
        return result;
    }

    std::vector<bool> min_cut(Int src) const {
        assert(0 <= src && src < Int(graph.size()));
        std::vector<bool> result(graph.size());
        result[src] = true;
        std::vector<Int> queue{src};
        for (std::size_t head = 0; head < queue.size(); ++head) {
            Int v = queue[head];
            for (auto e : graph[v])
                if (e.cap > Cap(0) && !result[e.dst]) {
                    result[e.dst] = true;
                    queue.push_back(e.dst);
                }
        }
        return result;
    }
};
}
}
