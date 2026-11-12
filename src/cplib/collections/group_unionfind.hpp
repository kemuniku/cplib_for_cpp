#pragma once
#include <cplib/common.hpp>

namespace cplib {
class UnionFind {
    std::vector<std::int32_t> par_or_siz;
    std::vector<Int> next, edge_cnt;

public:
    Int count;

    // グループの循環リストと辺数を持つ集合を構築する。O(n)。
    explicit UnionFind(Int n) : par_or_siz(n, -1), next(n), edge_cnt(n), count(n) {
        std::iota(next.begin(), next.end(), 0);
    }

    // 構築時の要素数を返す。連結成分数countとは異なる。O(1)。
    Int size() const {
        return par_or_siz.size();
    }

    // xの代表元を返して経路圧縮する。償却O(α(n))。
    Int root(Int x) {
        return par_or_siz[x] < 0 ? x : par_or_siz[x] = root(par_or_siz[x]);
    }

    // xとyが同じ集合か返す。償却O(α(n))。
    bool issame(Int x, Int y) {
        return root(x) == root(y);
    }

    // xとyの間に辺を追加し、その集合を結合する。償却O(α(n))。
    void unite(Int x, Int y) {
        x = root(x);
        y = root(y);
        if (x != y) {
            std::swap(next[x], next[y]);
            if (par_or_siz[x] > par_or_siz[y])
                std::swap(x, y);
            par_or_siz[x] += par_or_siz[y];
            par_or_siz[y] = x;
            edge_cnt[x] += edge_cnt[y];
            --count;
        }
        ++edge_cnt[x];
    }

    // xの集合の要素数を返す。償却O(α(n))。
    Int siz(Int x) {
        return -par_or_siz[root(x)];
    }

    // O(N)かけて、rootになっている頂点を列挙します。
    // 注意:O(root数)でないことに注意してください。
    std::vector<Int> roots() const {
        std::vector<Int> r;
        for (Int i = 0; i < Int(par_or_siz.size()); ++i)
            if (par_or_siz[i] < 0)
                r.push_back(i);
        return r;
    }

    // 指定頂点のグループを列挙する。O(グループの大きさ)。
    std::vector<Int> get_group(Int x) const {
        std::vector<Int> r;
        Int now = x;
        do {
            r.push_back(now);
            now = next[now];
        } while (now != x);
        return r;
    }

    // 全ての集合の要素を列挙する。O(n)。
    auto groups() const {
        std::vector<std::vector<Int>> r;
        for (Int x : roots())
            r.push_back(get_group(x));
        return r;
    }

    // 計算量は償却O(α(n))。
    // xの属するグループの辺の数を返します。
    Int edge_count(Int x) {
        return edge_cnt[root(x)];
    }

    // 計算量は償却O(α(n))。
    // xの属する連結成分が木かどうかを返します。
    bool is_tree(Int x) {
        return edge_count(x) == siz(x) - 1;
    }

    // 計算量は償却O(α(n))。
    // xの属する連結成分がなもりグラフかどうかを返します。
    bool is_namori(Int x) {
        return edge_count(x) == siz(x);
    }

    // 計算量は償却O(α(n))。
    // xの属する連結成分にサイクルがあるかどうかを返します。
    bool has_cycle(Int x) {
        return edge_count(x) >= siz(x);
    }

    // この構造の独立したコピーを返す。O(n)。
    UnionFind copy() const {
        return *this;
    }
};

// n要素のグループ管理Union-Findを生成する。O(n)。
inline UnionFind initUnionFind(Int n) {
    return UnionFind(n);
}

// xのグループの要素を列挙する。要素数をkとしてO(k)。
inline auto get_group(const UnionFind &u, Int x) {
    return u.get_group(x);
}

// 全てのグループの要素を列挙する。O(n)。
inline auto groups(const UnionFind &u) {
    return u.groups();
}

// xのグループの辺数を返す。償却O(α(n))。
inline Int edge_count(UnionFind &u, Int x) {
    return u.edge_count(x);
}

// xのグループが木か返す。償却O(α(n))。
inline bool is_tree(UnionFind &u, Int x) {
    return u.is_tree(x);
}

// xのグループがなもりグラフか返す。償却O(α(n))。
inline bool is_namori(UnionFind &u, Int x) {
    return u.is_namori(x);
}

// xのグループに閉路があるか返す。償却O(α(n))。
inline bool has_cycle(UnionFind &u, Int x) {
    return u.has_cycle(x);
}

// uの独立したコピーを返す。O(n)。
inline UnionFind copy(const UnionFind &u) {
    return u.copy();
}

// uのrootを呼び出す。計算量はu.rootに依存する。
template <class U, class... Args> auto root(U &u, Args... args) -> decltype(u.root(args...)) {
    return u.root(args...);
}

// uのissameを呼び出す。計算量はu.issameに依存する。
template <class U, class... Args> auto issame(U &u, Args... args) -> decltype(u.issame(args...)) {
    return u.issame(args...);
}

// uのuniteを呼び出す。計算量はu.uniteに依存する。
template <class U, class... Args> auto unite(U &u, Args... args) -> decltype(u.unite(args...)) {
    return u.unite(args...);
}

// uのsizを呼び出す。計算量はu.sizに依存する。
template <class U> auto siz(U &u, Int x) -> decltype(u.siz(x)) {
    return u.siz(x);
}

// uのrootsを呼び出す。計算量はu.rootsに依存する。
template <class U> auto roots(const U &u) -> decltype(u.roots()) {
    return u.roots();
}

}

namespace cplib {
using GroupUnionFind = cplib::UnionFind;
}
