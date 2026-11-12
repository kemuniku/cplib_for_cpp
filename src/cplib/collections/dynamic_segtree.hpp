#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <memory>

namespace cplib {
// 各座標を一つのノードに保持する動的セグメント木。Q回更新後の空間O(Q)。
template <class T> class DynamicSegmentTree {
    struct Node {
        Int index;
        T value, product;
        std::unique_ptr<Node> left, right;

        Node(Int i, T v) : index(i), value(v), product(v) {
        }
    };

    using Ptr = std::unique_ptr<Node>;
    Ptr root;
    Int length = 0, nodes = 0;
    std::function<T(T, T)> merge;
    T identity;

    static Ptr clone(const Ptr &p) {
        if (!p)
            return nullptr;
        auto q = std::make_unique<Node>(p->index, p->value);
        q->product = p->product;
        q->left = clone(p->left);
        q->right = clone(p->right);
        return q;
    }

    T product(const Ptr &node) const {
        return node ? node->product : identity;
    }

    void updateNode(Ptr &node, Int l, Int r, Int index, T value) {
        if (!node) {
            node = std::make_unique<Node>(index, value);
            ++nodes;
            return;
        }
        if (node->index == index)
            node->value = value;
        else {
            Int mid = l + (r - l) / 2;
            if (index < mid) {
                if (node->index < index) {
                    std::swap(node->index, index);
                    std::swap(node->value, value);
                }
                updateNode(node->left, l, mid, index, value);
            } else {
                if (index < node->index) {
                    std::swap(node->index, index);
                    std::swap(node->value, value);
                }
                updateNode(node->right, mid, r, index, value);
            }
        }
        node->product = merge(merge(product(node->left), node->value), product(node->right));
    }

    T getNode(const Ptr &node, Int l, Int r, Int ql, Int qr) const {
        if (!node || qr <= l || r <= ql)
            return identity;
        if (ql <= l && r <= qr)
            return node->product;
        Int mid = l + (r - l) / 2;
        T result = getNode(node->left, l, mid, ql, qr);
        if (ql <= node->index && node->index < qr)
            result = merge(result, node->value);
        return merge(result, getNode(node->right, mid, r, ql, qr));
    }

public:
    // 長さnの単位元列と演算を指定して空の木を構築する。O(1)。
    template <class Op> DynamicSegmentTree(Int n, Op op, T e) : length(n), merge(op), identity(e) {
        assert(n >= 0);
    }

    // 全ノードを複製する。ノード数をMとしてO(M)。
    DynamicSegmentTree(const DynamicSegmentTree &p)
        : root(clone(p.root)), length(p.length), nodes(p.nodes), merge(p.merge),
          identity(p.identity) {
    }

    // 木をムーブ構築する。O(1)。
    DynamicSegmentTree(DynamicSegmentTree &&) = default;
    // 木をムーブ代入する。代入先のノード数をMとしてO(M)。
    DynamicSegmentTree &operator=(DynamicSegmentTree &&) = default;

    // 木を複製して代入する。両者のノード数の合計をMとしてO(M)。
    DynamicSegmentTree &operator=(const DynamicSegmentTree &p) {
        if (this != &p) {
            DynamicSegmentTree copy(p);
            *this = std::move(copy);
        }
        return *this;
    }

    // 要素数を返す。O(1)。
    Int len() const {
        return length;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return len();
    }

    // 確保済みのノード数を返す。O(1)。
    Int node_count() const {
        return nodes;
    }

    // 一点更新O(log N)、新しく確保するノードは高々一つ。
    void update(Int index, T value) {
        assert(0 <= index && index < length);
        updateNode(root, 0, length, index, value);
    }

    // 半開区間[l,r)の積をO(log N)で返します。ノードは確保しません。
    T get(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= length);
        return l == r ? identity : getNode(root, 0, length, l, r);
    }

    // 両端を含む区間の積を返す。O(log(n+1))回の演算。
    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return get(resolve_index(length, s.a), resolve_index(length, s.b) + 1);
    }

    // 両端を含む区間の積を返す。O(log(n+1))回の演算。
    template <class L, class R> T operator[](ClosedSlice<L, R> s) const {
        return get(s);
    }

    // 1点の値をO(log N)で返します。未更新の座標では単位元を返します。
    T operator[](Int i) const {
        assert(0 <= i && i < length);
        auto node = root.get();
        while (node) {
            if (node->index == i)
                return node->value;
            node = i < node->index ? node->left.get() : node->right.get();
        }
        return identity;
    }

    struct Reference {
        DynamicSegmentTree *tree;
        Int index;

        // 対象の値を取得する。O(log(n+1))。
        operator T() const {
            return std::as_const(*tree)[index];
        }

        // 対象の値を設定する。O(log(n+1))回の演算。
        Reference &operator=(const T &v) {
            tree->update(index, v);
            return *this;
        }

        // 他の位置の値を読み取り代入する。O(log(n+1))回の演算。
        Reference &operator=(const Reference &r) {
            return *this = T(r);
        }
    };

    // i番目を読み書きする参照を返す。O(1)。
    Reference operator[](Int i) {
        return {this, i};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
    // 全区間の積を返す。O(1)。
    T get_all() const {
        return product(root);
    }
};

// 長さnの単位元列と演算を指定して空の木を生成する。O(1)。
template <class T, class Op> auto initDynamicSegmentTree(Int n, Op op, T e) {
    return DynamicSegmentTree<T>(n, op, e);
}

// 長さnの単位元列と演算を指定して空の木を生成する。O(1)。
template <class Op, class T> auto newDynamicSegWith(Int n, Op op, T e) {
    return initDynamicSegmentTree(n, op, e);
}

// 要素数を返す。O(1)。
template <class T> Int len(const DynamicSegmentTree<T> &s) {
    return s.len();
}

// 確保済みのノード数を返す。O(1)。
template <class T> Int node_count(const DynamicSegmentTree<T> &s) {
    return s.node_count();
}

// i番目の値を更新する。O(log(n+1))回の演算。
template <class T> void update(DynamicSegmentTree<T> &s, Int i, T value) {
    s.update(i, value);
}

// 半開区間[l,r)の積を返す。O(log(n+1))回の演算。
template <class T> T get(const DynamicSegmentTree<T> &s, Int l, Int r) {
    return s.get(l, r);
}

// 両端を含む区間の積を返す。O(log(n+1))回の演算。
template <class T, class L, class R> T get(const DynamicSegmentTree<T> &s, ClosedSlice<L, R> r) {
    return s.get(r);
}

// 全区間の積を返す。O(1)。
template <class T> T get_all(const DynamicSegmentTree<T> &s) {
    return s.get_all();
}
}
