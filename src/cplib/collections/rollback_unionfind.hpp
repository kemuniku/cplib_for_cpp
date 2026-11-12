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

    // 現在の集合数を返す。O(1)。
    Int count() const {
        return components;
    }

    // 結合操作の履歴数を状態として返す。O(1)。
    Int get_state() const {
        return history.size() / 2;
    }

    // xの代表元を返す。O(log n)。
    Int root(Int x) const {
        while (par_or_siz[x] >= 0)
            x = par_or_siz[x];
        return x;
    }

    // xとyが同じ集合に属するか返す。O(log n)。
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

    // 現在の状態を巻き戻し先として保存する。O(1)。
    void snapshot() {
        snap = get_state();
    }

    // 保存した巻き戻し先を初期状態に戻す。O(1)。
    void clear_snapshot() {
        snap = 0;
    }

    // 指定した状態まで戻す。省略時は保存した状態に戻す。取り消す操作数をkとしてO(k)。
    void rollback(Int state = -1) {
        if (state == -1)
            state = snap;
        assert(state >= 0 && state <= get_state());
        while (get_state() > state)
            undo();
    }

    // xの集合の要素数を返す。O(log n)。
    Int siz(Int x) const {
        return -par_or_siz[root(x)];
    }
};

// n要素の巻き戻し可能Union-Findを生成する。O(n)。
inline RollbackUnionFind initRollbackUnionFind(Int n) {
    return RollbackUnionFind(n);
}

// 現在の集合数を返す。O(1)。
inline Int count(const RollbackUnionFind &u) {
    return u.count();
}

// 現在の履歴数を状態として返す。O(1)。
inline Int get_state(const RollbackUnionFind &u) {
    return u.get_state();
}

// 直前の結合操作を取り消す。O(1)。
inline void undo(RollbackUnionFind &u) {
    u.undo();
}

// 現在の状態を巻き戻し先として保存する。O(1)。
inline void snapshot(RollbackUnionFind &u) {
    u.snapshot();
}

// 巻き戻し先を初期状態に戻す。O(1)。
inline void clear_snapshot(RollbackUnionFind &u) {
    u.clear_snapshot();
}

// 指定状態まで戻す。取り消す操作数をkとしてO(k)。
inline void rollback(RollbackUnionFind &u, Int state = -1) {
    u.rollback(state);
}
}
