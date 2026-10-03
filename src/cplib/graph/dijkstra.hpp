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

// 単一または複数始点から最短距離と前駆を返す。O((V+E) log(V+E))。
template <GraphTypes G, class Start, class T = typename G::cost_type>
RestoredDistances<T> restore_dijkstra(const G &g, const Start &start, T zero = T(0),
                                      T inf = detail::distance_inf<T>()) {
    std::priority_queue<std::pair<T, Int>, std::vector<std::pair<T, Int>>, std::greater<>> queue;
    RestoredDistances<T> result{std::vector<T>(g.len, inf), std::vector<Int>(g.len, -1)};
    auto init = [&](Int s) {
        queue.emplace(zero, s);
        result.costs[s] = zero;
    };
    if constexpr (std::is_integral_v<Start>)
        init(start);
    else
        for (Int s : start)
            init(s);
    while (!queue.empty()) {
        auto [cost, i] = queue.top();
        queue.pop();
        if (cost > result.costs[i])
            continue;
        for (auto [j, c] : g.to_and_cost(i)) {
            T temp = result.costs[i] + c;
            if (temp < result.costs[j]) {
                result.prev[j] = i;
                result.costs[j] = temp;
                queue.emplace(temp, j);
            }
        }
    }
    return result;
}

template <GraphTypes G, class Start, class T = typename G::cost_type>
auto dijkstra(const G &g, const Start &start, T zero = T(0), T inf = detail::distance_inf<T>()) {
    return restore_dijkstra(g, start, zero, inf).costs;
}

template <GraphTypes G, class T = typename G::cost_type>
ShortestPath<T> shortest_path_dijkstra(const G &g, Int start, Int goal, T zero = T(0),
                                       T inf = detail::distance_inf<T>()) {
    auto result = restore_dijkstra(g, start, zero, inf);
    return {restore_shortest_path_from_prev(result.prev, goal), result.costs[goal]};
}
}
