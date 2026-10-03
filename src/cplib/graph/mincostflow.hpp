#pragma once
#include <cplib/common.hpp>
#include <queue>
#include <stdexcept>
#include <type_traits>

namespace cplib {
template <class Cap, class Cost> struct MinCostFlowEdge {
    Int src, dst;
    Cap cap, flow;
    Cost cost;
};

template <class Cap = Int, class Cost = Int> class MinCostFlow {
    static_assert(std::is_integral_v<Cap> && std::is_integral_v<Cost> && std::is_signed_v<Cost>);

    struct Arc {
        Int dst, rev;
        Cap cap;
        Cost cost;
    };

    std::vector<std::vector<Arc>> graph;
    std::vector<std::pair<Int, Int>> positions;

public:
    // n頂点の最小費用流グラフを構築する。容量・費用型の省略時はInt。O(n)。
    explicit MinCostFlow(Int n) : graph(n) {
        assert(n >= 0);
    }

    // 容量cap、単位費用costの有向辺を追加し、辺番号を返す。償却O(1)。
    Int add_edge(Int src, Int dst, Cap cap, Cost cost) {
        assert(0 <= src && src < Int(graph.size()) && 0 <= dst && dst < Int(graph.size()) &&
               cap >= Cap(0) && cost != std::numeric_limits<Cost>::min());
        Int id = positions.size(), index = graph[src].size(),
            rev = graph[dst].size() + (src == dst);
        positions.emplace_back(src, index);
        graph[src].push_back({dst, rev, cap, cost});
        graph[dst].push_back({src, index, Cap(0), -cost});
        return id;
    }

    // i番目の辺の容量、現在の流量、単位費用を返す。O(1)。
    MinCostFlowEdge<Cap, Cost> get_edge(Int i) const {
        auto [src, index] = positions[i];
        auto e = graph[src][index];
        Cap f = graph[e.dst][e.rev].cap;
        return {src, e.dst, e.cap + f, f, e.cost};
    }

    // 追加順に全辺の情報を返す。O(E)。
    std::vector<MinCostFlowEdge<Cap, Cost>> get_edges() const {
        std::vector<MinCostFlowEdge<Cap, Cost>> out;
        out.reserve(positions.size());
        for (Int i = 0; i < Int(positions.size()); ++i)
            out.push_back(get_edge(i));
        return out;
    }

    // Bellman-Fordで初期ポテンシャルを求め、Dijkstraで増加路を選ぶ。
    // O(VE + A E log V)。到達可能な負閉路は例外。中間費用はCostに収まること。
    // 追加流量と最小費用の折れ点を返す。同じ傾きはまとめる。Aは増加回数。
    // 負費用辺に対応する。始点から到達可能な負閉路はstd::invalid_argument。
    std::vector<std::pair<Cap, Cost>> slope(Int src, Int dst,
                                            Cap limit = std::numeric_limits<Cap>::max()) {
        Int n = graph.size();
        assert(0 <= src && src < n && 0 <= dst && dst < n && src != dst && limit >= Cap(0));
        std::vector<std::pair<Cap, Cost>> out = {{0, 0}};
        if (limit == Cap(0))
            return out;
        std::vector<Cost> potential(n), distance(n);
        std::vector<bool> reached(n);
        std::vector<Int> prevVertex(n), prevEdge(n);
        reached[src] = true;
        for (Int phase = 0; phase < n; ++phase) {
            bool changed = false;
            for (Int v = 0; v < n; ++v)
                if (reached[v])
                    for (auto e : graph[v])
                        if (e.cap > Cap(0) &&
                            (!reached[e.dst] || potential[e.dst] > potential[v] + e.cost)) {
                            potential[e.dst] = potential[v] + e.cost;
                            reached[e.dst] = true;
                            changed = true;
                        }
            if (!changed)
                break;
            if (phase == n - 1)
                throw std::invalid_argument("始点から到達可能な負閉路があります");
        }
        Cap totalFlow = 0;
        Cost totalCost = 0, previousCost = 0;
        bool hasPrevious = false;
        while (totalFlow < limit) {
            std::fill(reached.begin(), reached.end(), false);
            distance[src] = 0;
            reached[src] = true;
            using Item = std::pair<Cost, Int>;
            std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
            queue.emplace(0, src);
            while (!queue.empty()) {
                auto [d, v] = queue.top();
                queue.pop();
                if (d != distance[v])
                    continue;
                for (Int i = 0; i < Int(graph[v].size()); ++i) {
                    auto e = graph[v][i];
                    if (e.cap == Cap(0))
                        continue;
                    Cost nd = d + (e.cost + potential[v] - potential[e.dst]);
                    if (!reached[e.dst] || nd < distance[e.dst]) {
                        reached[e.dst] = true;
                        distance[e.dst] = nd;
                        prevVertex[e.dst] = v;
                        prevEdge[e.dst] = i;
                        queue.emplace(nd, e.dst);
                    }
                }
            }
            if (!reached[dst])
                break;
            for (Int v = 0; v < n; ++v)
                if (reached[v])
                    potential[v] += distance[v];
            Cost unitCost = potential[dst] - potential[src];
            Cap pushed = limit - totalFlow;
            for (Int v = dst; v != src; v = prevVertex[v])
                pushed = std::min(pushed, graph[prevVertex[v]][prevEdge[v]].cap);
            Cost nextCost = totalCost;
            if (unitCost != Cost(0))
                nextCost += Cost(pushed) * unitCost;
            for (Int v = dst; v != src; v = prevVertex[v]) {
                Int u = prevVertex[v], i = prevEdge[v], rev = graph[u][i].rev;
                graph[u][i].cap -= pushed;
                graph[v][rev].cap += pushed;
            }
            totalFlow += pushed;
            totalCost = nextCost;
            if (hasPrevious && previousCost == unitCost)
                out.pop_back();
            out.emplace_back(totalFlow, totalCost);
            previousCost = unitCost;
            hasPrevious = true;
        }
        return out;
    }

    // limit以下の流量を追加し、追加流量と費用を返す。計算量と制約はslopeと同じ。
    std::pair<Cap, Cost> flow(Int src, Int dst, Cap limit = std::numeric_limits<Cap>::max()) {
        return slope(src, dst, limit).back();
    }
};

template <class Cap = Int, class Cost = Int>
auto initMinCostFlow(Int n, Cap = Cap(0), Cost = Cost(0)) {
    return MinCostFlow<Cap, Cost>(n);
}

template <class Cap, class Cost>
Int add_edge(MinCostFlow<Cap, Cost> &g, Int s, Int t, Cap cap, Cost cost) {
    return g.add_edge(s, t, cap, cost);
}

template <class Cap, class Cost> auto get_edge(const MinCostFlow<Cap, Cost> &g, Int i) {
    return g.get_edge(i);
}

template <class Cap, class Cost> auto get_edges(const MinCostFlow<Cap, Cost> &g) {
    return g.get_edges();
}

template <class Cap, class Cost>
auto slope(MinCostFlow<Cap, Cost> &g, Int s, Int t, Cap limit = std::numeric_limits<Cap>::max()) {
    return g.slope(s, t, limit);
}

template <class Cap, class Cost>
auto flow(MinCostFlow<Cap, Cost> &g, Int s, Int t, Cap limit = std::numeric_limits<Cap>::max()) {
    return g.flow(s, t, limit);
}
}
