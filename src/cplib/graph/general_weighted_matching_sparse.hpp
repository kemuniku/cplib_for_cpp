#pragma once
#include <cplib/graph/internal/weighted_matching_engine.hpp>

namespace cplib {
// Positive edges only; parallel edges use the largest weight.
// 最大重みと、それを達成する頂点ペア列(u < v)を返す。時間O(V+E+KM log K)、追加領域O(V+E)。
// Kは正の非ループ辺に接する頂点数、Mは多重辺をまとめた後の正の非ループ辺数。
// 辺数の最大化は保証しない。自己ループと重み0以下の辺は無視し、多重辺は最大重みを使う。
// 重みは符号付き整数型で、総和はInt。K>0のとき最大正重みWはstd::numeric_limits<Int>::max() / (4*K)以下とする。
// 入力は変更しない。静的グラフはbuild()が必要。
template <UnDirectedGraph G> auto maximum_weight_matching_sparse(const G &g) {
    return detail::weighted_matching::solve<true>(g);
}
}
