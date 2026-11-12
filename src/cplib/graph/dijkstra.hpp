#pragma once
#include <cplib/graph/private/warshall_floyd_common.hpp>
#include <cplib/graph/restore_shortest_path_from_prev.hpp>
#include <queue>

namespace cplib {
template <class T> struct RestoredDistances {
    std::vector<T> costs;
    std::vector<Int> prev;
};

template <class T> struct ShortestPath {
    std::vector<Int> path;
    T cost;
};

// 更新関数addと負・0・正を返す比較cmpで距離と親を返す。距離更新は単調非減少であること。
// 始点数Sとして時間 O(V+(E+S)log(E+S+2))、領域 O(V+E+S)。演算・比較はO(1)とする。
template <GraphTypes G, class Start, class Add, class T, class Compare>
RestoredDistances<T> restore_dijkstra(const G &g, const Start &start, Add add, T zero, T inf, Compare cmp) {
    using Entry = std::pair<T, Int>;
    auto lowerPriority = [cmp](const Entry &a, const Entry &b) {
        auto order = cmp(a.first, b.first);
        return order == 0 ? a.second > b.second : order > 0;
    };
    std::priority_queue<Entry, std::vector<Entry>, decltype(lowerPriority)> queue(lowerPriority);
    RestoredDistances<T> result{std::vector<T>(g.size(), inf), std::vector<Int>(g.size(), -1)};
    std::vector<bool> initialized(g.size());
    auto init = [&](Int vertex) {
        if (!initialized[vertex]) {
            initialized[vertex] = true;
            result.costs[vertex] = zero;
            queue.emplace(zero, vertex);
        }
    };
    if constexpr (std::is_integral_v<Start>)
        init(start);
    else
        for (Int vertex : start)
            init(vertex);
    while (!queue.empty()) {
        auto [cost, vertex] = queue.top();
        queue.pop();
        if (cmp(result.costs[vertex], cost) < 0)
            continue;
        for (auto [next, weight] : g.to_and_cost(vertex)) {
            T candidate = add(result.costs[vertex], weight);
            if (cmp(candidate, result.costs[next]) < 0) {
                result.prev[next] = vertex;
                result.costs[next] = candidate;
                queue.emplace(candidate, next);
            }
        }
    }
    return result;
}
// 更新関数addと通常の比較で距離と親を返す。時間 O(V+(E+S)log(E+S+2))、領域 O(V+E+S)。
template <GraphTypes G, class Start, class Add, class T>
auto restore_dijkstra(const G &g, const Start &start, Add add, T zero, T inf) {
    return restore_dijkstra(g, start, add, zero, inf, [](const T &a, const T &b) {
        return a < b ? -1 : b < a ? 1 : 0;
    });
}
// 非負辺の通常の加算で距離と親を返す。時間 O(V+(E+S)log(E+S+2))、領域 O(V+E+S)。
template <GraphTypes G, class Start, class T = typename G::cost_type>
auto restore_dijkstra(const G &g, const Start &start, T zero = T(0), T inf = detail::distance_inf<T>()) {
    return restore_dijkstra(g, start, [](const T &distance, const auto &edge) { return distance + edge; }, zero, inf);
}
// 更新関数addと比較cmpで各頂点への距離を返す。時間 O(V+(E+S)log(E+S+2))、領域 O(V+E+S)。
template <GraphTypes G, class Start, class Add, class T, class Compare>
auto dijkstra(const G &g, const Start &start, Add add, T zero, T inf, Compare cmp) {
    return restore_dijkstra(g, start, add, zero, inf, cmp).costs;
}
// 更新関数addと通常の比較で距離を返す。時間 O(V+(E+S)log(E+S+2))、領域 O(V+E+S)。
template <GraphTypes G, class Start, class Add, class T>
auto dijkstra(const G &g, const Start &start, Add add, T zero, T inf) {
    return restore_dijkstra(g, start, add, zero, inf).costs;
}
// addとcmpで経路と距離を返す。到達不能なら経路は終点のみ。時間 O(V+E log(E+2))、領域 O(V+E)。
template <GraphTypes G, class Add, class T, class Compare>
ShortestPath<T> shortest_path_dijkstra(const G &g, Int start, Int goal, Add add, T zero, T inf, Compare cmp) {
    auto result = restore_dijkstra(g, start, add, zero, inf, cmp);
    return {restore_shortest_path_from_prev(result.prev, goal), result.costs[goal]};
}
// addと通常の比較で経路と距離を返す。時間 O(V+E log(E+2))、領域 O(V+E)。
template <GraphTypes G, class Add, class T>
ShortestPath<T> shortest_path_dijkstra(const G &g, Int start, Int goal, Add add, T zero, T inf) {
    auto result = restore_dijkstra(g, start, add, zero, inf);
    return {restore_shortest_path_from_prev(result.prev, goal), result.costs[goal]};
}

// 単一または複数始点からの最短距離を返す。非負重みのDijkstra法。始点数をSとしてO(V+(S+1)E+(E+S) log(E+S+1))。
template <GraphTypes G, class Start, class T = typename G::cost_type>
auto dijkstra(const G &g, const Start &start, T zero = T(0), T inf = detail::distance_inf<T>()) {
    return restore_dijkstra(g, start, zero, inf).costs;
}

// 始点から終点への経路と距離を返す。到達不能時の経路は終点のみ。非負重みのDijkstra法でO(V+E log(E+2))。
template <GraphTypes G, class T = typename G::cost_type>
ShortestPath<T> shortest_path_dijkstra(const G &g, Int start, Int goal, T zero = T(0),
                                       T inf = detail::distance_inf<T>()) {
    auto result = restore_dijkstra(g, start, zero, inf);
    return {restore_shortest_path_from_prev(result.prev, goal), result.costs[goal]};
}
}
