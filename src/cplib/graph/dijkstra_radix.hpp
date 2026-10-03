#pragma once
#include <cplib/collections/radix_heap.hpp>
#include <cplib/graph/dijkstra.hpp>

namespace cplib {
namespace detail {
template <class T> constexpr T radix_distance_inf() {
    if constexpr (std::is_same_v<T, Int>)
        return T(INF64);
    else if constexpr (std::is_same_v<T, std::int32_t>)
        return T(INF32);
    else
        return std::numeric_limits<T>::max();
}

template <GraphTypes G, class Start, class T>
RestoredDistances<T> radix_dijkstra(const G &g, const Start &start, T zero, T inf, Int goal = -1) {
    static_assert(std::is_integral_v<T>);
    RadixHeap<T, Int> q(zero);
    RestoredDistances<T> out{std::vector<T>(g.len, inf), std::vector<Int>(g.len, -1)};
    auto init = [&](Int s) {
        assert(0 <= s && s < g.len);
        if (out.costs[s] != zero) {
            out.costs[s] = zero;
            q.push(zero, s);
        }
    };
    if constexpr (std::is_integral_v<Start>)
        init(start);
    else
        for (Int s : start)
            init(s);
    while (q.len()) {
        auto [cost, u] = q.pop();
        if (cost > out.costs[u])
            continue;
        if (u == goal)
            break;
        for (auto [v, w] : g.to_and_cost(u)) {
            assert(w >= 0);
            if (cost > std::numeric_limits<T>::max() - w)
                continue;
            T d = cost + w;
            if (d < out.costs[v]) {
                out.costs[v] = d;
                out.prev[v] = u;
                q.push(d, v);
            }
        }
    }
    return out;
}
}

// 非負整数重みの最短距離と直前の頂点を返す。同距離の経路は不定。
// 重みなしグラフでは各辺を1とする。Sは始点数、Bは距離型のビット数としてO(V + (E + S)B)。
template <GraphTypes G, class Start, class T = typename G::cost_type>
auto restore_dijkstra_radix(const G &g, const Start &s, T zero = T(0),
                            T inf = detail::radix_distance_inf<T>()) {
    return detail::radix_dijkstra(g, s, zero, inf);
}

// 非負整数重みの最短距離を返す。inf以上は未到達扱い。重みなしグラフでは各辺を1とする。
// 始点sは一頂点または頂点列。Sは始点数、Bは距離型のビット数としてO(V + (E + S)B)。
template <GraphTypes G, class Start, class T = typename G::cost_type>
auto dijkstra_radix(const G &g, const Start &s, T zero = T(0),
                    T inf = detail::radix_distance_inf<T>()) {
    return restore_dijkstra_radix(g, s, zero, inf).costs;
}

// 目的地tが確定するまで探索する。未到達時はpath={t}、cost=inf。
// 非負整数重みを指定し、重みなしグラフでは各辺を1とする。Bは距離型のビット数としてO(V + EB)。
template <GraphTypes G, class T = typename G::cost_type>
ShortestPath<T> shortest_path_dijkstra_radix(const G &g, Int s, Int t, T zero = T(0),
                                             T inf = detail::radix_distance_inf<T>()) {
    auto r = detail::radix_dijkstra(g, s, zero, inf, t);
    return {restore_shortest_path_from_prev(r.prev, t), r.costs[t]};
}
}
