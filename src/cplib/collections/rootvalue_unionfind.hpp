#pragma once
#include <cplib/collections/unionfind.hpp>
#include <functional>

namespace cplib {
template <class T> class RootValueUnionFind {
    std::vector<std::int32_t> par_or_siz;
    std::function<void(T &, T &)> op;

public:
    Int count;
    std::vector<T> values;

    // 初期値と結合時の操作を指定してn要素を構築する。O(n)。
    RootValueUnionFind(Int n, std::function<void(T &, T &)> operation, std::vector<T> initial)
        : par_or_siz(n, -1), op(std::move(operation)), count(n), values(std::move(initial)) {
    }

    // xの代表元を返して経路圧縮する。償却O(α(n))。
    Int root(Int x) {
        return par_or_siz[x] < 0 ? x : par_or_siz[x] = root(par_or_siz[x]);
    }

    // xとyが同じ集合に属するか返す。償却O(α(n))。
    bool issame(Int x, Int y) {
        return root(x) == root(y);
    }

    // 大きい側の値をopで更新して結合する。償却O(α(n))+opの計算量。
    void unite(Int x, Int y) {
        x = root(x);
        y = root(y);
        if (x == y)
            return;
        if (par_or_siz[x] > par_or_siz[y])
            std::swap(x, y);
        op(values[x], values[y]);
        par_or_siz[x] += par_or_siz[y];
        par_or_siz[y] = x;
        --count;
    }

    // xの集合の要素数を返す。償却O(α(n))。
    Int siz(Int x) {
        return -par_or_siz[root(x)];
    }

    // xの集合に保持された値への参照を返す。償却O(α(n))。
    T &get(Int x) {
        return values[root(x)];
    }

    // xの集合の値を設定する。償却O(α(n))と値の代入に要する時間。
    void set(Int x, T value) {
        values[root(x)] = std::move(value);
    }
};

// 全要素に同じ初期値を与えて構築する。値のコピーがO(1)ならO(n)。
template <class T, class Op> auto initRootValueUnionFind(Int n, Op op, const T &value) {
    return RootValueUnionFind<T>(n, op, std::vector<T>(n, value));
}

// 配列から初期値を与えて構築する。値のコピーがO(1)ならO(n)。
template <class T, class Op> auto initRootValueUnionFind(Int n, Op op, std::span<const T> values) {
    return RootValueUnionFind<T>(n, op, std::vector<T>(values.begin(), values.end()));
}

// 配列から初期値を与えて構築する。値のコピーがO(1)ならO(n)。
template <class T, class Op>
auto initRootValueUnionFind(Int n, Op op, const std::vector<T> &values) {
    return RootValueUnionFind<T>(n, op, values);
}

// factoryで各初期値を生成して構築する。O(n)とfactoryのn回分の計算量。
template <class T, class Op, class Factory>
    requires std::invocable<Factory>
auto initRootValueUnionFind(Int n, Op op, Factory factory) {
    std::vector<T> values;
    values.reserve(n);
    for (Int i = 0; i < n; ++i)
        values.push_back(factory());
    return RootValueUnionFind<T>(n, op, std::move(values));
}

// xの集合の値への参照を返す。償却O(α(n))。
template <class T> T &get(RootValueUnionFind<T> &u, Int x) {
    return u.get(x);
}

// xの集合の値を設定する。償却O(α(n))と値の代入に要する時間。
template <class T> void setValue(RootValueUnionFind<T> &u, Int x, T value) {
    u.set(x, std::move(value));
}
}
