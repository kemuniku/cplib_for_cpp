#pragma once
#include <cplib/collections/unionfind.hpp>

namespace cplib {
template <class T = Int> class WeightedUnionFind {
    std::vector<std::int32_t> par_or_siz;
    std::vector<T> potential_diff;

public:
    Int count;

    // ポテンシャル付き集合を構築する。O(n)。
    explicit WeightedUnionFind(Int n) : par_or_siz(n, -1), potential_diff(n, T(0)), count(n) {
    }

    // 構築時の要素数を返す。連結成分数countとは異なる。O(1)。
    Int size() const {
        return par_or_siz.size();
    }

    // xの代表元を返して経路圧縮し、ポテンシャル差を更新する。償却O(α(n))。
    Int root(Int x) {
        if (par_or_siz[x] < 0)
            return x;
        Int p = par_or_siz[x], r = root(p);
        potential_diff[x] += potential_diff[p];
        par_or_siz[x] = r;
        return r;
    }

    // xの代表元に対するポテンシャルを返す。償却O(α(n))。
    T potential(Int x) {
        root(x);
        return potential_diff[x];
    }

    // xとyが同じ集合に属するか返す。償却O(α(n))。
    bool issame(Int x, Int y) {
        return root(x) == root(y);
    }

    // 計算量は償却O(α(n))。
    // 同じ連結成分内のpotential[y]-potential[x]を返します
    T diff(Int x, Int y) {
        Int rx = root(x), ry = root(y);
        assert(rx == ry);
        (void)rx;
        (void)ry;
        return potential_diff[y] - potential_diff[x];
    }

    // potential[y]-potential[x]=wを満たすよう結合する。矛盾ならfalse。償却O(α(n))。
    bool unite(Int x, Int y, T w) {
        w = w + potential(x) - potential(y);
        x = root(x);
        y = root(y);
        if (x == y)
            return w == T(0);
        if (par_or_siz[x] > par_or_siz[y]) {
            std::swap(x, y);
            w = -w;
        }
        par_or_siz[x] += par_or_siz[y];
        par_or_siz[y] = x;
        --count;
        potential_diff[y] = w;
        return true;
    }

    // xの集合の要素数を返す。償却O(α(n))。
    Int siz(Int x) {
        return -par_or_siz[root(x)];
    }
};

// n要素の重み付きUnion-Findを生成する。O(n)。
template <class T = Int> auto initWeightedUnionFind(Int n) {
    return WeightedUnionFind<T>(n);
}

// xの代表元に対するポテンシャルを返す。償却O(α(n))。
template <class T> T potential(WeightedUnionFind<T> &u, Int x) {
    return u.potential(x);
}

// 同じ集合のyとxのポテンシャル差を返す。償却O(α(n))。
template <class T> T diff(WeightedUnionFind<T> &u, Int x, Int y) {
    return u.diff(x, y);
}
}
