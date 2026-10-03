#pragma once
#include <cplib/common.hpp>

namespace cplib {
class UnionFind {
    std::vector<std::int32_t> par_or_siz;

public:
    Int count;

    // n個の集合を構築する。O(n)。
    explicit UnionFind(Int n) : par_or_siz(n, -1), count(n) {
    }

    // 経路圧縮して根を返す。償却O(α(n))。
    Int root(Int x) {
        return par_or_siz[x] < 0 ? x : par_or_siz[x] = root(par_or_siz[x]);
    }

    bool issame(Int x, Int y) {
        return root(x) == root(y);
    }

    void unite(Int x, Int y) {
        x = root(x);
        y = root(y);
        if (x == y)
            return;
        if (par_or_siz[x] > par_or_siz[y])
            std::swap(x, y);
        par_or_siz[x] += par_or_siz[y];
        par_or_siz[y] = x;
        --count;
    }

    Int siz(Int x) {
        return -par_or_siz[root(x)];
    }

    // 根を列挙する。O(n)。
    std::vector<Int> roots() const {
        std::vector<Int> r;
        r.reserve(count);
        for (Int i = 0; i < Int(par_or_siz.size()); ++i)
            if (par_or_siz[i] < 0)
                r.push_back(i);
        return r;
    }

    UnionFind copy() const {
        return *this;
    }
};

inline UnionFind initUnionFind(Int n) {
    return UnionFind(n);
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

inline UnionFind copy(const UnionFind &u) {
    return u.copy();
}
}
