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

    // 構築時の要素数を返す。連結成分数countとは異なる。O(1)。
    Int size() const {
        return par_or_siz.size();
    }

    // xの代表元を返し、経路圧縮する。償却O(α(n))。
    Int root(Int x) {
        return par_or_siz[x] < 0 ? x : par_or_siz[x] = root(par_or_siz[x]);
    }

    // xとyが同じ集合に属するか返す。償却O(α(n))。
    bool issame(Int x, Int y) {
        return root(x) == root(y);
    }

    // xとyの集合を結合する。償却O(α(n))。
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

    // xの集合の要素数を返す。償却O(α(n))。
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

    // この構造の独立したコピーを返す。O(n)。
    UnionFind copy() const {
        return *this;
    }
};

// n個の要素を持つUnion-Findを生成する。O(n)。
inline UnionFind initUnionFind(Int n) {
    return UnionFind(n);
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

// uの独立したコピーを返す。O(n)。
inline UnionFind copy(const UnionFind &u) {
    return u.copy();
}
}
