#pragma once
#include <cplib/graph/private/warshall_floyd_common.hpp>
#include <cplib/graph/private/warshall_floyd_avx_kernel.hpp>
namespace cplib::detail {
struct Avx2Warshall {
template<bool Check,class T> static bool run(std::vector<std::vector<T>>& d,T zero,T inf){
if constexpr(std::is_same_v<T,Int> || std::is_same_v<T,std::int32_t>){
std::vector<T*> rows;rows.reserve(d.size());for(auto& row:d)rows.push_back(row.data());
if constexpr(std::is_same_v<T,Int>){
if constexpr(Check)return avx2_kernel::cplib_warshall_floyd_int64_avx2(rows.data(),d.size(),zero,inf);
else { avx2_kernel::cplib_warshall_floyd_nonnegative_int64_avx2(rows.data(),d.size(),zero,inf);return false;}
}
else{
if constexpr(Check)return avx2_kernel::cplib_warshall_floyd_int32_avx2(rows.data(),d.size(),zero,inf);
else { avx2_kernel::cplib_warshall_floyd_nonnegative_int32_avx2(rows.data(),d.size(),zero,inf);return false;}
}
} else return ScalarWarshall::run<Check>(d,zero,inf);
}
};
}
namespace cplib {
// 全点対最短路を上書きする。O(V^3)。
template<class T> void warshall_floyd_inplace(std::vector<std::vector<T>>& d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){cplib::detail::warshall_inplace<cplib::detail::Avx2Warshall,true>(d,zero,inf);}
// 行列から距離を求める。O(V^3)時間・O(V^2)領域。
template<class T> auto warshall_floyd(std::vector<std::vector<T>> d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){ cplib::warshall_floyd_inplace(d,zero,inf);return d;}
// グラフから距離を求める。O(V^3)時間・O(V^2)領域。
template<cplib::GraphTypes G,class T=typename G::cost_type> auto warshall_floyd(const G& g,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){auto d=cplib::detail::graph_distance_matrix(g,zero,inf); cplib::warshall_floyd_inplace(d,zero,inf);return d;}
// 全点対最短路を上書きする。O(V^3)。
template<class T> void warshall_floyd_nonnegative_inplace(std::vector<std::vector<T>>& d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){cplib::detail::warshall_inplace<cplib::detail::Avx2Warshall,false>(d,zero,inf);}
// 行列から距離を求める。O(V^3)時間・O(V^2)領域。
template<class T> auto warshall_floyd_nonnegative(std::vector<std::vector<T>> d,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){ cplib::warshall_floyd_nonnegative_inplace(d,zero,inf);return d;}
// グラフから距離を求める。O(V^3)時間・O(V^2)領域。
template<cplib::GraphTypes G,class T=typename G::cost_type> auto warshall_floyd_nonnegative(const G& g,T zero=T(0),T inf=cplib::detail::distance_inf<T>()){auto d=cplib::detail::graph_distance_matrix(g,zero,inf); cplib::warshall_floyd_nonnegative_inplace(d,zero,inf);return d;}
}
