#pragma once
#include <cplib/collections/rootvalue_unionfind.hpp>

namespace cplib {
class DynamicBipartite {
    Int N;
    RootValueUnionFind<Int> uf;
    bool bipartite = true;

    static auto make_uf(Int n) {
        std::vector<Int> values(2 * n);
        std::fill(values.begin() + n, values.end(), 1);
        return initRootValueUnionFind<Int>(2 * n, [](Int &x, Int &y) { x += y; }, values);
    }

public:
    Int cnt_sum;

    explicit DynamicBipartite(Int n) : N(n), uf(make_uf(n)), cnt_sum(n) {
    }

    // 色の表裏を2頂点で保持する。各操作償却O(α(N))。
    // 辺(u,v)を新たに追加したときにグラフが二部グラフかどうか判定する
    // 既に二部グラフではない場合は常にfalseが返る
    bool can_unite(Int u, Int v) {
        return bipartite && !uf.issame(u, v);
    }

    // 辺(u,v)の追加
    void unite(Int u, Int v) {
        if (uf.issame(u, v + N) || !bipartite)
            return;
        if (uf.issame(u, v)) {
            bipartite = false;
            cnt_sum = -1;
            return;
        }
        u = uf.root(u);
        v = uf.root(v);
        cnt_sum -= std::max(uf.get(u), uf.get(u + N)) + std::max(uf.get(v), uf.get(v + N));
        uf.unite(u, N + v);
        uf.unite(v, N + u);
        cnt_sum += std::max(uf.get(u), uf.get(u + N));
    }

    bool is_bipartite() const {
        return bipartite;
    }

    bool issame(Int u, Int v) {
        return uf.issame(u, v) || uf.issame(u, v + N);
    }
};

inline auto initDynamicBipartite(Int n) {
    return DynamicBipartite(n);
}

inline bool can_unite(DynamicBipartite &g, Int u, Int v) {
    return g.can_unite(u, v);
}

inline bool is_bipartite(const DynamicBipartite &g) {
    return g.is_bipartite();
}
}
