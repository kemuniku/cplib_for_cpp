#pragma once
#include <cplib/graph/private/warshall_floyd_common.hpp>
namespace cplib {
// 全点対最短路を上書きする。O(V^3)。
template<class T> void warshall_floyd_inplace(std::vector<std::vector<T>>& d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){cplib::detail::warshall_inplace<cplib::detail::ScalarWarshall,true>(d,zero,inf);}
// 行列から距離を求める。O(V^3)時間・O(V^2)領域。
template<class T> auto warshall_floyd(std::vector<std::vector<T>> d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){ cplib::warshall_floyd_inplace(d,zero,inf);return d;}
// グラフから距離を求める。O(V^3)時間・O(V^2)領域。
template<cplib::GraphTypes G,class T=typename G::cost_type> auto warshall_floyd(const G& g,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){auto d=cplib::detail::graph_distance_matrix(g,zero,inf); cplib::warshall_floyd_inplace(d,zero,inf);return d;}
// 全点対最短路を上書きする。O(V^3)。
template<class T> void warshall_floyd_nonnegative_inplace(std::vector<std::vector<T>>& d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){cplib::detail::warshall_inplace<cplib::detail::ScalarWarshall,false>(d,zero,inf);}
// 行列から距離を求める。O(V^3)時間・O(V^2)領域。
template<class T> auto warshall_floyd_nonnegative(std::vector<std::vector<T>> d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){ cplib::warshall_floyd_nonnegative_inplace(d,zero,inf);return d;}
// グラフから距離を求める。O(V^3)時間・O(V^2)領域。
template<cplib::GraphTypes G,class T=typename G::cost_type> auto warshall_floyd_nonnegative(const G& g,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){auto d=cplib::detail::graph_distance_matrix(g,zero,inf); cplib::warshall_floyd_nonnegative_inplace(d,zero,inf);return d;}
}
