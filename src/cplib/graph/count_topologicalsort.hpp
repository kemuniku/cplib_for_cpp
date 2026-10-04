#pragma once
#include <cplib/graph/topologicalsort.hpp>

namespace cplib {
// トポロジカル順序の個数を返す。時間O(V * 2^V + E)、空間O(2^V + V)。
// 空グラフは1、閉路を含む場合は0。多重辺は同じ制約として扱い、重みは無視する。
// 静的グラフは事前にbuildが必要。頂点数は20程度までを想定する。
// V < 63かつ答えがIntに収まることが必要（V <= 20なら収まる）。
template <DirectedGraph G> Int count_topologicalsort(const G &g) {
    Int n = g.len;
    assert(n < 63);
    if (!isDAG(g))
        return 0;
    std::vector<UInt> pred(n);
    for (Int u = 0; u < n; ++u)
        for (auto [v, c] : g.to_and_cost(u))
            pred[v] |= UInt(1) << u;
    UInt size = UInt(1) << n;
    std::vector<Int> dp(size);
    dp[0] = 1;
    for (UInt mask = 0; mask < size; ++mask)
        if (dp[mask])
            for (Int v = 0; v < n; ++v) {
                UInt bit = UInt(1) << v;
                if (!(mask & bit) && (mask & pred[v]) == pred[v])
                    dp[mask | bit] += dp[mask];
            }
    return dp.back();
}
}
