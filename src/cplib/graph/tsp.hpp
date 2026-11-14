#pragma once
#include <cplib/graph/maxk_dijkstra.hpp>

namespace cplib {
namespace detail {
template <class T>
auto solve_tsp(std::vector<std::vector<T>> dist, const std::vector<T> &start, T inf, bool floyd) {
    Int n = dist.size();
    if (floyd)
        for (Int k = 0; k < n; ++k)
            for (Int i = 0; i < n; ++i)
                for (Int j = 0; j < n; ++j)
                    dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
    std::vector<std::vector<T>> dp(n, std::vector<T>(Int(1) << n, inf));
    for (Int i = 0; i < n; ++i)
        dp[i][Int(1) << i] = start[i];
    for (Int bit = 0; bit < (Int(1) << n); ++bit)
        for (Int i = 0; i < n; ++i)
            if ((bit & (Int(1) << i)) && dp[i][bit] != inf)
                for (Int j = 0; j < n; ++j)
                    if (!(bit & (Int(1) << j)))
                        dp[j][bit | (Int(1) << j)] =
                            std::min(dp[j][bit | (Int(1) << j)], dp[i][bit] + dist[i][j]);
    return std::pair{std::move(dp), std::move(dist)};
}
}

// 全頂点訪問の部分集合DP。O(N²2^N)、floyd指定時は先に距離閉包を取る。
// start_vからすべての頂点を訪問して、終点はどこでもよいときのコストの最小値。
template <class T>
T tspPathCostFrom(const std::vector<std::vector<T>> &dist, Int start_v, T zero, T inf,
                  bool floydwarshall = true) {
    std::vector<T> start(dist.size(), inf);
    start[start_v] = zero;
    auto [dp, d] = detail::solve_tsp(dist, start, inf, floydwarshall);
    T output_value = inf;
    for (const auto &row : dp)
        output_value = std::min(output_value, row.back());
    return output_value;
}

// 指定した始点から全頂点を訪問する最小費用を返す。O(N²2^N + N³)。
template <class T>
T tspPathCostFrom(const std::vector<std::vector<T>> &dist, Int start_v, bool floydwarshall = true) {
    return tspPathCostFrom(dist, start_v, T(0), detail::distance_inf<T>(), floydwarshall);
}

// 計算量O(N²2^N + N³)。
// 元実装どおり、最後の訪問点からgoalへ移動する費用も加える（goalの再訪可）。
// start_vからすべての頂点を訪問して、終点をgoal_vに指定したときのコストの最小値。
template <class T>
T tspPathCostFromTo(const std::vector<std::vector<T>> &dist, Int start_v, Int goal_v, T zero, T inf,
                    bool floydwarshall = true) {
    std::vector<T> start(dist.size(), inf);
    start[start_v] = zero;
    auto [dp, d] = detail::solve_tsp(dist, start, inf, floydwarshall);
    T output_value = inf;
    for (Int i = 0; i < Int(dist.size()); ++i)
        output_value = std::min(output_value, dp[i].back() + d[i][goal_v]);
    return output_value;
}

// 指定した始点から全頂点を訪問して終点に至る最小費用を返す。O(N²2^N + N³)。
template <class T>
T tspPathCostFromTo(const std::vector<std::vector<T>> &dist, Int start_v, Int goal_v,
                    bool floydwarshall = true) {
    return tspPathCostFromTo(dist, start_v, goal_v, T(0), detail::distance_inf<T>(), floydwarshall);
}

// 計算量O(N²2^N + N³)。
// グラフのどこかからスタートし、すべての頂点を通るまでのコストの最小値。
template <class T>
T tspPathAnyStart(const std::vector<std::vector<T>> &dist, T zero, T inf,
                  bool floydwarshall = true) {
    std::vector<T> start(dist.size(), zero);
    auto [dp, d] = detail::solve_tsp(dist, start, inf, floydwarshall);
    T output_value = inf;
    for (const auto &row : dp)
        output_value = std::min(output_value, row.back());
    return output_value;
}

// 任意の始点から全頂点を訪問する最小費用を返す。O(N²2^N + N³)。
template <class T>
T tspPathAnyStart(const std::vector<std::vector<T>> &dist, bool floydwarshall = true) {
    return tspPathAnyStart(dist, T(0), detail::distance_inf<T>(), floydwarshall);
}

// 指定頂点数KとしてO(K(V+E)log(V+E)+K²)。辺重みは非負とする。
// 指定頂点間の最短距離を辺重みにする。重みなし版は循環バケットBFSを使用。
template <WeightedGraph G, class T = typename G::cost_type>
auto toContractionGraph(const G &g, std::span<const Int> vertices, T zero = T(0),
                        T inf = detail::distance_inf<T>()) {
    auto output_value = initWeightedDirectedGraph<T>(vertices.size());
    for (Int i = 0; i < Int(vertices.size()); ++i) {
        auto res = dijkstra(g, vertices[i], zero, inf);
        for (Int j = 0; j < Int(vertices.size()); ++j)
            output_value.add_edge(i, j, res[vertices[j]]);
    }
    return output_value;
}

// 指定頂点間の最短距離を重みにする。指定頂点数KとしてO(K(V+E)+K²)。
template <UnWeightedGraph G> auto toContractionGraph(const G &g, std::span<const Int> vertices) {
    auto output_value = initWeightedDirectedGraph<Int>(vertices.size());
    for (Int i = 0; i < Int(vertices.size()); ++i) {
        auto res = maxk_dijkstra(g, vertices[i], Int(1));
        for (Int j = 0; j < Int(vertices.size()); ++j)
            output_value.add_edge(i, j, res[vertices[j]]);
    }
    return output_value;
}

// グラフの隣接行列を返す。O(V²+E)。
// 自己辺がなければ対角もnoneのまま。多重辺は最小費用を採用。
template <DirectedGraph G, class T = typename G::cost_type>
    requires WeightedGraph<G>
auto to_adjacency_matrix(const G &g, T none = detail::distance_inf<T>()) {
    std::vector<std::vector<T>> output_value(g.len, std::vector<T>(g.len, none));
    for (Int i = 0; i < g.len; ++i)
        for (auto [j, c] : g.to_and_cost(i))
            output_value[i][j] = std::min(output_value[i][j], T(c));
    return output_value;
}
}
