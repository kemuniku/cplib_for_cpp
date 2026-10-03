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

    Int root(Int x) {
        return par_or_siz[x] < 0 ? x : par_or_siz[x] = root(par_or_siz[x]);
    }

    bool issame(Int x, Int y) {
        return root(x) == root(y);
    }

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

    auto groups() const {
        std::vector<std::vector<Int>> r;
        for (Int x : roots())
            r.push_back(get_group(x));
        return r;
    }

    // xの属するグループの辺の数を返します。
    Int edge_count(Int x) {
        return edge_cnt[root(x)];
    }

    // xの属する連結成分が木かどうかを返します。
    bool is_tree(Int x) {
        return edge_count(x) == siz(x) - 1;
    }

    // xの属する連結成分がなもりグラフかどうかを返します。
    bool is_namori(Int x) {
        return edge_count(x) == siz(x);
    }

    // xの属する連結成分にサイクルがあるかどうかを返します。
    bool has_cycle(Int x) {
        return edge_count(x) >= siz(x);
    }

    UnionFind copy() const {
        return *this;
    }
};

inline UnionFind initUnionFind(Int n) {
    return UnionFind(n);
}

inline auto get_group(const UnionFind &u, Int x) {
    return u.get_group(x);
}

inline auto groups(const UnionFind &u) {
    return u.groups();
}

inline Int edge_count(UnionFind &u, Int x) {
    return u.edge_count(x);
}

inline bool is_tree(UnionFind &u, Int x) {
    return u.is_tree(x);
}

inline bool is_namori(UnionFind &u, Int x) {
    return u.is_namori(x);
}

inline bool has_cycle(UnionFind &u, Int x) {
    return u.has_cycle(x);
}

inline UnionFind copy(const UnionFind &u) {
    return u.copy();
}

template <class U, class... Args> auto root(U &u, Args... args) -> decltype(u.root(args...)) {
    return u.root(args...);
}

template <class U, class... Args> auto issame(U &u, Args... args) -> decltype(u.issame(args...)) {
    return u.issame(args...);
}

template <class U, class... Args> auto unite(U &u, Args... args) -> decltype(u.unite(args...)) {
    return u.unite(args...);
}

template <class U> auto siz(U &u, Int x) -> decltype(u.siz(x)) {
    return u.siz(x);
}

template <class U> auto roots(const U &u) -> decltype(u.roots()) {
    return u.roots();
}

}

namespace cplib {
using GroupUnionFind = cplib::UnionFind;
}
