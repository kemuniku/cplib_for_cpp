#pragma once
#include <cplib/common.hpp>
#include <queue>
#include <functional>

namespace cplib {
struct MinCostBFlowEdge {
    Int src, dst, lower, upper, cost, flow = 0;
};
struct MinCostBFlowResult {
    bool feasible = false;
    __int128 cost = 0;
};
class MinCostBFlow {
    struct Arc {
        Int dst;
        __int128 cap, cost;
    };
    std::vector<__int128> supply_, potential_;
    std::vector<MinCostBFlowEdge> edges_;
    bool solved_ = false;
public:
    // n頂点を作る。正のbは供給（流出−流入=b）。時間・領域 O(n)。
    explicit MinCostBFlow(Int n) : supply_(n) {
        assert(n >= 0);
    }
    // 頂点数を返す。計算量 O(1)。
    Int size() const {
        return supply_.size();
    }
    // 供給量を加算する。負値は需要。計算量 O(1)。
    void add_supply(Int v, Int amount) {
        assert(0 <= v && v < size());
        supply_[v] += amount;
        solved_ = false;
    }
    // 上下限・費用を持つ辺を追加する。負の上下限や自己ループも許す。償却 O(1)。
    Int add_edge(Int src, Int dst, Int lower, Int upper, Int cost) {
        assert(0 <= src && src < size() && 0 <= dst && dst < size() && lower <= upper);
        Int id = edges_.size();
        edges_.push_back({src, dst, lower, upper, cost});
        solved_ = false;
        return id;
    }
    // 容量スケーリングで毎回解き直す。実行不能ならfeasible=false、cost=0。
    // 時間 O((N+M)² log(N+M+2) log(U+2)+NM)、領域 O(N+M)。
    // Uは上下限差・下限適用後の収支の絶対値の最大値。中間値は符号付き128ビットに収まること。
    MinCostBFlowResult solve() {
        solved_ = false;
        Int n = size();
        auto balance = supply_;
        __int128 total = 0, bound = 1;
        for (auto b : balance)
            total += b;
        if (total != 0)
            return {};
        std::vector<std::vector<Int>> graph(n);
        std::vector<Arc> arcs;
        for (const auto &e : edges_) {
            __int128 width = __int128(e.upper) - e.lower;
            graph[e.src].push_back(arcs.size());
            arcs.push_back({e.dst, width, e.cost});
            graph[e.dst].push_back(arcs.size());
            arcs.push_back({e.src, 0, -__int128(e.cost)});
            balance[e.src] -= e.lower;
            balance[e.dst] += e.lower;
            bound = std::max(bound, width);
        }
        for (auto b : balance)
            bound = std::max(bound, b < 0 ? -b : b);
        __int128 delta = 1;
        while (delta <= bound / 2)
            delta *= 2;
        std::vector<__int128> potential(n), distance(n);
        std::vector<bool> reached(n);
        std::vector<Int> parent(n);
        for (; delta > 0; delta /= 2) {
            for (Int u = 0; u < n; ++u) {
                for (Int i : graph[u]) {
                    Int v = arcs[i].dst;
                    if (arcs[i].cost + potential[u] - potential[v] < 0) {
                        __int128 amount = arcs[i].cap - arcs[i].cap % delta;
                        arcs[i].cap -= amount;
                        arcs[i ^ 1].cap += amount;
                        balance[u] -= amount;
                        balance[v] += amount;
                    }
                }
            }
            while (true) {
                using Entry = std::pair<__int128, Int>;
                std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
                for (Int v = 0; v < n; ++v) {
                    reached[v] = balance[v] >= delta;
                    parent[v] = -1;
                    if (reached[v]) {
                        distance[v] = 0;
                        queue.emplace(0, v);
                    }
                }
                Int target = -1;
                while (!queue.empty()) {
                    auto [d, u] = queue.top();
                    queue.pop();
                    if (d != distance[u])
                        continue;
                    if (balance[u] <= -delta) {
                        target = u;
                        break;
                    }
                    for (Int i : graph[u]) {
                        if (arcs[i].cap < delta)
                            continue;
                        Int v = arcs[i].dst;
                        __int128 nd = d + arcs[i].cost + potential[u] - potential[v];
                        if (!reached[v] || nd < distance[v]) {
                            reached[v] = true;
                            distance[v] = nd;
                            parent[v] = i;
                            queue.emplace(nd, v);
                        }
                    }
                }
                if (target < 0)
                    break;
                for (Int v = 0; v < n; ++v)
                    potential[v] += reached[v] ? std::min(distance[v], distance[target]) : distance[target];
                __int128 amount = -balance[target];
                Int root = target;
                while (parent[root] >= 0) {
                    Int i = parent[root];
                    amount = std::min(amount, arcs[i].cap);
                    root = arcs[i ^ 1].dst;
                }
                amount = std::min(amount, balance[root]);
                amount -= amount % delta;
                for (Int v = target; v != root;) {
                    Int i = parent[v];
                    arcs[i].cap -= amount;
                    arcs[i ^ 1].cap += amount;
                    v = arcs[i ^ 1].dst;
                }
                balance[root] -= amount;
                balance[target] += amount;
            }
        }
        for (auto b : balance)
            if (b != 0)
                return {};
        potential_.assign(n, 0);
        for (Int phase = 0; phase < n; ++phase) {
            bool changed = false;
            for (Int u = 0; u < n; ++u) {
                for (Int i : graph[u]) {
                    Int v = arcs[i].dst;
                    if (arcs[i].cap > 0 && potential_[v] > potential_[u] + arcs[i].cost) {
                        potential_[v] = potential_[u] + arcs[i].cost;
                        changed = true;
                    }
                }
            }
            if (!changed)
                break;
        }
        MinCostBFlowResult result{true, 0};
        for (Int i = 0; i < Int(edges_.size()); ++i) {
            __int128 flow = __int128(edges_[i].lower) + arcs[2 * i + 1].cap;
            edges_[i].flow = Int(flow);
            result.cost += flow * edges_[i].cost;
        }
        solved_ = true;
        return result;
    }
    // 解が有効なとき元辺と流量を返す。計算量 O(1)。
    MinCostBFlowEdge get_edge(Int i) const {
        assert(solved_ && 0 <= i && i < Int(edges_.size()));
        return edges_[i];
    }
    // 解が有効なとき全辺と流量をコピーして返す。時間・領域 O(M)。
    auto get_edges() const {
        assert(solved_);
        return edges_;
    }
    // 残余辺の被約費用が非負となる双対解を返す。時間・領域 O(N)。
    auto get_potential() const {
        assert(solved_);
        return potential_;
    }
};
// n頂点の最小費用b-flowを作る。時間・領域 O(n)。
inline auto initMinCostBFlow(Int n) {
    return MinCostBFlow(n);
}
}
