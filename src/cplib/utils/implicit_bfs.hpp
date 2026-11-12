#pragma once
#include <cplib/utils/constants.hpp>
#include <functional>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace cplib {
template <class T> using ImplicitBfsAdjacent = std::function<std::vector<T>(T)>;
template <class T> using ImplicitBfsFinish = std::function<bool(T)>;
template <class T> struct ImplicitBfsResult {
    std::unordered_map<T, Int> costs;
    std::unordered_map<T, T> prev;
};
template <class T> struct ImplicitBfsPath {
    std::vector<T> path;
    Int cost;
};
namespace detail {
template <class T> inline constexpr bool isBfsStartList = false;
template <class T, class Allocator>
inline constexpr bool isBfsStartList<std::vector<T, Allocator>> = true;

template <class T, class Adjacent, class Finish>
auto implicitBfsRun(const std::vector<T> &start, Adjacent adjacent, Finish finish, Int inf) {
    ImplicitBfsResult<T> result;
    ImplicitBfsPath<T> path{{}, inf};
    std::queue<T> queue;
    for (const T &s : start)
        if (result.costs.emplace(s, 0).second)
            queue.push(s);
    while (!queue.empty()) {
        T vertex = queue.front();
        queue.pop();
        Int cost = result.costs.at(vertex);
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
        for (const T &next : adjacent(vertex)) {
            if (result.costs.emplace(next, cost + 1).second) {
                result.prev.emplace(next, vertex);
                queue.push(next);
            }
        }
    }
    return std::pair{std::move(result), std::move(path)};
}
}
// 重複を除いた各始点からBFSし距離と親を返す。期待時間 O(S+V+E)、領域 O(V)。
template <class T, class Adjacent>
auto restore_implicit_bfs(const std::vector<T> &start, Adjacent adjacent) {
    return detail::implicitBfsRun(start, adjacent, [](const T &) { return false; }, INF64).first;
}
// 単一始点の距離と親を返す。期待時間 O(V+E)、領域 O(V)。
template <class T, class Adjacent>
auto restore_implicit_bfs(const T &start, Adjacent adjacent) {
    return restore_implicit_bfs(std::vector<T>{start}, adjacent);
}
// 各始点からの最小辺数を返す。期待時間 O(S+V+E)、領域 O(V)。
template <class T, class Adjacent>
auto implicit_bfs(const std::vector<T> &start, Adjacent adjacent) {
    return restore_implicit_bfs(start, adjacent).costs;
}
// 単一始点からの最小辺数を返す。期待時間 O(V+E)、領域 O(V)。
template <class T, class Adjacent>
auto implicit_bfs(const T &start, Adjacent adjacent) {
    return implicit_bfs(std::vector<T>{start}, adjacent);
}
// 初めてfinishを満たす頂点までの距離を返す。不在はinf。期待時間 O(S+V+E)、領域 O(V)。
template <class T, class Adjacent, class Finish>
    requires (std::is_invocable_v<Adjacent, const T &> && std::is_invocable_r_v<bool, Finish, T>)
Int implicit_bfs_until(const std::vector<T> &start, Adjacent adjacent, Finish finish, Int inf = INF64) {
    return detail::implicitBfsRun(start, adjacent, finish, inf).second.cost;
}
// 単一始点の途中終了BFS。不在はinf。期待時間 O(V+E)、領域 O(V)。
template <class T, class Adjacent, class Finish>
    requires (!detail::isBfsStartList<T> &&
              std::is_invocable_v<Adjacent, const T &> && std::is_invocable_r_v<bool, Finish, T>)
Int implicit_bfs_until(const T &start, Adjacent adjacent, Finish finish, Int inf = INF64) {
    return implicit_bfs_until(std::vector<T>{start}, adjacent, finish, inf);
}
// いずれかの終点への最短経路を返す。不在は空経路とinf。期待時間 O(S+G+V+E)、領域 O(V+G)。
template <class T, class Adjacent>
auto shortest_path_implicit_bfs(const std::vector<T> &start, const std::vector<T> &goal,
                                Adjacent adjacent, Int inf = INF64) {
    if (goal.empty())
        return ImplicitBfsPath<T>{{}, inf};
    std::unordered_set<T> goals(goal.begin(), goal.end());
    return detail::implicitBfsRun(start, adjacent, [&](const T &v) { return goals.contains(v); }, inf).second;
}
// 多始点・単一終点の最短経路を返す。期待時間 O(S+V+E)、領域 O(V)。
template <class T, class Adjacent>
auto shortest_path_implicit_bfs(const std::vector<T> &start, const T &goal,
                                Adjacent adjacent, Int inf = INF64) {
    return shortest_path_implicit_bfs(start, std::vector<T>{goal}, adjacent, inf);
}
// 単一始点・多終点の最短経路を返す。期待時間 O(G+V+E)、領域 O(V+G)。
template <class T, class Adjacent>
auto shortest_path_implicit_bfs(const T &start, const std::vector<T> &goal,
                                Adjacent adjacent, Int inf = INF64) {
    return shortest_path_implicit_bfs(std::vector<T>{start}, goal, adjacent, inf);
}
// 単一始点・単一終点の最短経路を返す。期待時間 O(V+E)、領域 O(V)。
template <class T, class Adjacent>
auto shortest_path_implicit_bfs(const T &start, const T &goal, Adjacent adjacent, Int inf = INF64) {
    return shortest_path_implicit_bfs(std::vector<T>{start}, std::vector<T>{goal}, adjacent, inf);
}
// いずれかの終点への最小距離を返す。期待時間 O(S+G+V+E)、領域 O(V+G)。
template <class T, class Adjacent>
Int implicit_bfs_until(const std::vector<T> &start, Adjacent adjacent,
                       const std::vector<T> &goal, Int inf = INF64) {
    return shortest_path_implicit_bfs(start, goal, adjacent, inf).cost;
}
// 単一始点からいずれかの終点への最小距離を返す。期待時間 O(G+V+E)、領域 O(V+G)。
template <class T, class Adjacent>
Int implicit_bfs_until(const T &start, Adjacent adjacent, const std::vector<T> &goal, Int inf = INF64) {
    return implicit_bfs_until(std::vector<T>{start}, adjacent, goal, inf);
}
}
