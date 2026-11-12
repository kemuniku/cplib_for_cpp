#pragma once
#include <cplib/utils/constants.hpp>
#include <functional>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace cplib {
template <class T> using ImplicitDijkstraAdjacent = std::function<std::vector<std::pair<T, Int>>(T)>;
template <class T> using ImplicitDijkstraFinish = std::function<bool(T)>;
template <class T> struct ImplicitDijkstraResult {
    std::unordered_map<T, Int> costs;
    std::unordered_map<T, T> prev;
};
template <class T> struct ImplicitDijkstraPath {
    std::vector<T> path;
    Int cost;
};
namespace detail {
template <class T> inline constexpr bool isDijkstraStartList = false;
template <class T, class Allocator>
inline constexpr bool isDijkstraStartList<std::vector<T, Allocator>> = true;

template <class T, class Adjacent, class Finish>
auto implicitDijkstraRun(const std::vector<T> &start, Adjacent adjacent, Finish finish, Int inf) {
    ImplicitDijkstraResult<T> result;
    ImplicitDijkstraPath<T> path{{}, inf};
    std::unordered_map<T, Int> ids;
    std::vector<T> vertices;
    using Entry = std::pair<Int, Int>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    for (const T &s : start) {
        if (result.costs.emplace(s, 0).second) {
            Int id = vertices.size();
            ids.emplace(s, id);
            vertices.push_back(s);
            queue.emplace(0, id);
        }
    }
    while (!queue.empty()) {
        auto [cost, id] = queue.top();
        queue.pop();
        T vertex = vertices[id];
        if (cost != result.costs.at(vertex))
            continue;
        if (finish(vertex)) {
            path.cost = cost;
            path.path.push_back(vertex);
            auto it = result.prev.find(vertex);
            while (it != result.prev.end()) {
                vertex = it->second;
                path.path.push_back(vertex);
                it = result.prev.find(vertex);
            }
            std::reverse(path.path.begin(), path.path.end());
            break;
        }
        for (const auto &[next, weight] : adjacent(vertex)) {
            assert(weight >= 0);
            Int nextCost = cost + weight;
            auto found = result.costs.find(next);
            if (found == result.costs.end() || nextCost < found->second) {
                if (!ids.contains(next)) {
                    ids.emplace(next, vertices.size());
                    vertices.push_back(next);
                }
                result.costs.insert_or_assign(next, nextCost);
                result.prev.insert_or_assign(next, vertex);
                queue.emplace(nextCost, ids.at(next));
            }
        }
    }
    return std::pair{std::move(result), std::move(path)};
}
}
// 重複を除いた各始点から非負辺のDijkstra探索し距離と親を返す。期待時間 O((S+V+E) log(S+V+E+2))、領域 O(V+E)。
template <class T, class Adjacent>
auto restore_implicit_dijkstra(const std::vector<T> &start, Adjacent adjacent) {
    return detail::implicitDijkstraRun(start, adjacent, [](const T &) { return false; }, INF64).first;
}
// 単一始点の距離と親を返す。期待時間 O((V+E) log(V+E+2))、領域 O(V+E)。
template <class T, class Adjacent>
auto restore_implicit_dijkstra(const T &start, Adjacent adjacent) {
    return restore_implicit_dijkstra(std::vector<T>{start}, adjacent);
}
// 各始点からの最小費用を返す。期待時間 O((S+V+E) log(S+V+E+2))、領域 O(V+E)。
template <class T, class Adjacent>
auto implicit_dijkstra(const std::vector<T> &start, Adjacent adjacent) {
    return restore_implicit_dijkstra(start, adjacent).costs;
}
// 単一始点からの最小費用を返す。期待時間 O((V+E) log(V+E+2))、領域 O(V+E)。
template <class T, class Adjacent>
auto implicit_dijkstra(const T &start, Adjacent adjacent) {
    return implicit_dijkstra(std::vector<T>{start}, adjacent);
}
// 初めてfinishを満たす頂点までの距離を返す。不在はinf。期待時間 O((S+V+E) log(S+V+E+2))、領域 O(V+E)。
template <class T, class Adjacent, class Finish>
    requires (std::is_invocable_v<Adjacent, const T &> && std::is_invocable_r_v<bool, Finish, T>)
Int implicit_dijkstra_until(const std::vector<T> &start, Adjacent adjacent, Finish finish, Int inf = INF64) {
    return detail::implicitDijkstraRun(start, adjacent, finish, inf).second.cost;
}
// 単一始点の途中終了非負辺のDijkstra探索。不在はinf。期待時間 O((V+E) log(V+E+2))、領域 O(V+E)。
template <class T, class Adjacent, class Finish>
    requires (!detail::isDijkstraStartList<T> &&
              std::is_invocable_v<Adjacent, const T &> && std::is_invocable_r_v<bool, Finish, T>)
Int implicit_dijkstra_until(const T &start, Adjacent adjacent, Finish finish, Int inf = INF64) {
    return implicit_dijkstra_until(std::vector<T>{start}, adjacent, finish, inf);
}
// いずれかの終点への最短経路を返す。不在は空経路とinf。期待時間 O(G+(S+V+E) log(S+V+E+2))、領域 O(V+E+G)。
template <class T, class Adjacent>
auto shortest_path_implicit_dijkstra(const std::vector<T> &start, const std::vector<T> &goal,
                                Adjacent adjacent, Int inf = INF64) {
    if (goal.empty())
        return ImplicitDijkstraPath<T>{{}, inf};
    std::unordered_set<T> goals(goal.begin(), goal.end());
    return detail::implicitDijkstraRun(start, adjacent, [&](const T &v) { return goals.contains(v); }, inf).second;
}
// 多始点・単一終点の最短経路を返す。期待時間 O((S+V+E) log(S+V+E+2))、領域 O(V+E)。
template <class T, class Adjacent>
auto shortest_path_implicit_dijkstra(const std::vector<T> &start, const T &goal,
                                Adjacent adjacent, Int inf = INF64) {
    return shortest_path_implicit_dijkstra(start, std::vector<T>{goal}, adjacent, inf);
}
// 単一始点・多終点の最短経路を返す。期待時間 O(G+(V+E) log(V+E+2))、領域 O(V+E+G)。
template <class T, class Adjacent>
auto shortest_path_implicit_dijkstra(const T &start, const std::vector<T> &goal,
                                Adjacent adjacent, Int inf = INF64) {
    return shortest_path_implicit_dijkstra(std::vector<T>{start}, goal, adjacent, inf);
}
// 単一始点・単一終点の最短経路を返す。期待時間 O((V+E) log(V+E+2))、領域 O(V+E)。
template <class T, class Adjacent>
auto shortest_path_implicit_dijkstra(const T &start, const T &goal, Adjacent adjacent, Int inf = INF64) {
    return shortest_path_implicit_dijkstra(std::vector<T>{start}, std::vector<T>{goal}, adjacent, inf);
}
// いずれかの終点への最小距離を返す。期待時間 O(G+(S+V+E) log(S+V+E+2))、領域 O(V+E+G)。
template <class T, class Adjacent>
Int implicit_dijkstra_until(const std::vector<T> &start, Adjacent adjacent,
                       const std::vector<T> &goal, Int inf = INF64) {
    return shortest_path_implicit_dijkstra(start, goal, adjacent, inf).cost;
}
// 単一始点からいずれかの終点への最小距離を返す。期待時間 O(G+(V+E) log(V+E+2))、領域 O(V+E+G)。
template <class T, class Adjacent>
Int implicit_dijkstra_until(const T &start, Adjacent adjacent, const std::vector<T> &goal, Int inf = INF64) {
    return implicit_dijkstra_until(std::vector<T>{start}, adjacent, goal, inf);
}
}
