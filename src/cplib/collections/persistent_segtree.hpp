#pragma once
#include <cplib/collections/persistent_lazysegtree.hpp>

namespace cplib {
template <class T = Int> struct SegmentTreeNode {
    T value;
    std::shared_ptr<SegmentTreeNode> left, right;
};

namespace detail {
struct PersistentNoAction {};
}

template <class T = Int> class PSegmentTree {
    PersistentLazySegmentTree<T, detail::PersistentNoAction> tree;

    explicit PSegmentTree(decltype(tree) t) : tree(std::move(t)) {
    }

public:
    // vから永続セグメント木をO(N)で構築する。opは結合演算、eは単位元。
    template <class V, class Op>
    PSegmentTree(const V &v, Op op, T e)
        : tree(
              v, op, e, [](auto, T x) { return x; },
              [](auto, auto) { return detail::PersistentNoAction{}; },
              detail::PersistentNoAction{}) {
    }

    // 要素数をO(1)で返す。
    Int len() const {
        return tree.len();
    }

    // pをvに置き換えた新しい版を、時間・追加領域O(log N)で返す。
    PSegmentTree update(Int p, T v) const {
        return PSegmentTree(tree.update(p, std::move(v)));
    }

    // 半開区間[l,r)をsrcの同じ区間で置き換えた新しい版を返す。木の操作は時間・追加領域O(log N)。
    // 両方の木は同じ長さ・演算・単位元で構築する。元の版は変更しない。
    // 独立に構築した木の部分コピーでは領域所有権も共有する。K個の木の統合に伴う領域参照の移動は全操作でO(K log K)。
    PSegmentTree copy_range(const PSegmentTree &src, Int l, Int r) const {
        return PSegmentTree(tree.copy_range(src.tree, l, r));
    }

    // 閉区間sをsrcの同じ区間で置き換えた新しい版を返す。
    template <class L, class R>
    PSegmentTree copy_range(const PSegmentTree &src, ClosedSlice<L, R> s) const {
        return PSegmentTree(tree.copy_range(src.tree, s));
    }

    // 半開区間[l,r)の積をO(log N)で返す。
    T get(Int l, Int r) const {
        return tree.get(l, r);
    }

    // getと同じく半開区間[l,r)の積をO(log N)で返す。
    T query(Int l, Int r) const {
        return get(l, r);
    }

    // 全要素の積をO(1)で返す。空なら単位元。
    T get_all() const {
        return tree.get_all();
    }

    // iの要素をO(log N)で返す。
    template <class I> T operator[](I i) const {
        return tree[i];
    }

    // 閉区間sの積をO(log N)で返す。
    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return tree.get(s);
    }

    // 変数を更新後の版へO(log N)で差し替える。他の変数に保存した版は変化しない。
    template <class I> void set(I i, T v) {
        tree.set(i, std::move(v));
    }

    // f(get(l,r))を満たす最大のrをO(log N)で返す。fは区間の拡大に対して単調で、単位元に対してtrueを返すこと。
    template <class Pred> Int max_right(Int l, Pred f) const {
        return tree.max_right(l, f);
    }

    // f(get(l,r))を満たす最小のlをO(log N)で返す。fは区間の拡大に対して単調で、単位元に対してtrueを返すこと。
    template <class Pred> Int min_left(Int r, Pred f) const {
        return tree.min_left(r, f);
    }

    // 要素を空白区切りで文字列化する。O(N + 出力長)。
    std::string str() const {
        return tree.str();
    }
};
template <class V, class Op, class T> PSegmentTree(const V &, Op, T) -> PSegmentTree<T>;
template <class T> using PersistentSegmentTree = PSegmentTree<T>;

template <class V, class Op, class T> auto initPersistentSegmentTree(const V &v, Op op, T e) {
    return PSegmentTree<T>(v, op, e);
}

template <class V, class Op, class T> auto initSegmentTree(const V &v, Op op, T e) {
    return initPersistentSegmentTree(v, op, e);
}

template <class T> auto update(const PSegmentTree<T> &s, Int i, T v) {
    return s.update(i, v);
}

template <class T> T query(const PSegmentTree<T> &s, Int l, Int r) {
    return s.query(l, r);
}

template <class... A> auto newPersistentSegWith(A &&...args) {
    return initPersistentSegmentTree(std::forward<A>(args)...);
}

template <class... A> auto newSegWith(A &&...args) {
    return initPersistentSegmentTree(std::forward<A>(args)...);
}
}
