#pragma once
#include <cplib/graph/private/warshall_floyd_common.hpp>

namespace cplib {
// 全点対最短路を上書きする。O(V^3)。
// 正方隣接行列を最短距離で上書きする。到達不能はinf、負閉路を経由できる組は-inf。O(V^3)。
template <class T>
void warshall_floyd_inplace(std::vector<std::vector<T>> &d, T zero = T(0),
                            T inf = cplib::detail::distance_inf<T>()) {
    cplib::detail::warshall_inplace<cplib::detail::ScalarWarshall, true>(d, zero, inf);
}

// 行列から距離を求める。O(V^3)時間・O(V^2)領域。
// 隣接行列から全点対最短路をO(V^3)で返す。到達不能はinf、負閉路を経由できる組は-inf。入力は変更しない。
template <class T>
auto warshall_floyd(std::vector<std::vector<T>> d, T zero = T(0),
                    T inf = cplib::detail::distance_inf<T>()) {
    cplib::warshall_floyd_inplace(d, zero, inf);
    return d;
}

// グラフから距離を求める。O(V^3)時間・O(V^2)領域。
// 全点対最短路を返す。到達不能はinf、負閉路を経由できる組は-inf。O(V^3)。
template <cplib::GraphTypes G, class T = typename G::cost_type>
auto warshall_floyd(const G &g, T zero = T(0), T inf = cplib::detail::distance_inf<T>()) {
    auto d = cplib::detail::graph_distance_matrix(g, zero, inf);
    cplib::warshall_floyd_inplace(d, zero, inf);
    return d;
}

// 全点対最短路を上書きする。O(V^3)。
// 負閉路がない隣接行列を最短距離で上書きする。辺なしはinf、対角はzero。負辺も許容し、負閉路検査は行わない。O(V^3)。
template <class T>
void warshall_floyd_nonnegative_inplace(std::vector<std::vector<T>> &d, T zero = T(0),
                                        T inf = cplib::detail::distance_inf<T>()) {
    cplib::detail::warshall_inplace<cplib::detail::ScalarWarshall, false>(d, zero, inf);
}

// 行列から距離を求める。O(V^3)時間・O(V^2)領域。
// 非負辺の正方隣接行列から距離行列を求める。時間O(V^3)、追加領域O(V^2)。入力は変更しない。
template <class T>
auto warshall_floyd_nonnegative(std::vector<std::vector<T>> d, T zero = T(0),
                                T inf = cplib::detail::distance_inf<T>()) {
    cplib::warshall_floyd_nonnegative_inplace(d, zero, inf);
    return d;
}

// グラフから距離を求める。O(V^3)時間・O(V^2)領域。
// 非負辺のグラフから距離行列を求める。時間O(V^3)、追加領域O(V^2)。
template <cplib::GraphTypes G, class T = typename G::cost_type>
auto warshall_floyd_nonnegative(const G &g, T zero = T(0),
                                T inf = cplib::detail::distance_inf<T>()) {
    auto d = cplib::detail::graph_distance_matrix(g, zero, inf);
    cplib::warshall_floyd_nonnegative_inplace(d, zero, inf);
    return d;
}
}
