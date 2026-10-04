#pragma once
#include <cplib/collections/unionfind.hpp>

namespace cplib {
class RollbackUnionFind {
    Int components, snap = 0;
    std::vector<std::int32_t> par_or_siz;
    std::vector<std::pair<Int, std::int32_t>> history;

public:
    // 巻き戻し可能な集合を初期化する。O(n)。
    explicit RollbackUnionFind(Int n) : components(n), par_or_siz(n, -1) {
    }

    // 構築時の要素数を返す。count()は現在の連結成分数。O(1)。
    Int size() const {
        return par_or_siz.size();
    }

    Int count() const {
        return components;
    }

    Int get_state() const {
        return history.size() / 2;
    }

    Int root(Int x) const {
        while (par_or_siz[x] >= 0)
            x = par_or_siz[x];
        return x;
    }

    bool issame(Int x, Int y) const {
        return root(x) == root(y);
    }

    // 結合の成否にかかわらず履歴を追加する。O(log n)。
    // 2頂点の集合を結合し、結合できたときだけ成分数を減らす。O(log N)。
    bool unite(Int x, Int y) {
        x = root(x);
        y = root(y);
        auto sx = par_or_siz[x], sy = par_or_siz[y];
        history.emplace_back(x, sx);
        history.emplace_back(y, sy);
        if (x == y)
            return false;
        if (sx > sy)
            std::swap(x, y);
        par_or_siz[x] += par_or_siz[y];
        par_or_siz[y] = x;
        --components;
        return true;
    }

    // 直前の結合操作を取り消し、成分数も復元する。O(1)。
    void undo() {
        assert(!history.empty());
        if (history.back().first != history[history.size() - 2].first)
            ++components;
        for (int i = 0; i < 2; ++i) {
            auto [x, s] = history.back();
            history.pop_back();
            par_or_siz[x] = s;
        }
        if (snap > get_state())
            snap = 0;
    }

    void snapshot() {
        snap = get_state();
    }

    void clear_snapshot() {
        snap = 0;
    }

    void rollback(Int state = -1) {
        if (state == -1)
            state = snap;
        assert(state >= 0 && state <= get_state());
        while (get_state() > state)
            undo();
    }

    Int siz(Int x) const {
        return -par_or_siz[root(x)];
    }
};

inline RollbackUnionFind initRollbackUnionFind(Int n) {
    return RollbackUnionFind(n);
}

inline Int count(const RollbackUnionFind &u) {
    return u.count();
}

inline Int get_state(const RollbackUnionFind &u) {
    return u.get_state();
}

inline void undo(RollbackUnionFind &u) {
    u.undo();
}

inline void snapshot(RollbackUnionFind &u) {
    u.snapshot();
}

inline void clear_snapshot(RollbackUnionFind &u) {
    u.clear_snapshot();
}

inline void rollback(RollbackUnionFind &u, Int state = -1) {
    u.rollback(state);
}
}
