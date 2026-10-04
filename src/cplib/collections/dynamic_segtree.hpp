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
    template <class Op> DynamicSegmentTree(Int n, Op op, T e) : length(n), merge(op), identity(e) {
        assert(n >= 0);
    }

    DynamicSegmentTree(const DynamicSegmentTree &p)
        : root(clone(p.root)), length(p.length), nodes(p.nodes), merge(p.merge),
          identity(p.identity) {
    }

    DynamicSegmentTree(DynamicSegmentTree &&) = default;
    DynamicSegmentTree &operator=(DynamicSegmentTree &&) = default;

    DynamicSegmentTree &operator=(const DynamicSegmentTree &p) {
        if (this != &p) {
            DynamicSegmentTree copy(p);
            *this = std::move(copy);
        }
        return *this;
    }

    Int len() const {
        return length;
    }

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

    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return get(resolve_index(length, s.a), resolve_index(length, s.b) + 1);
    }

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

        operator T() const {
            return std::as_const(*tree)[index];
        }

        Reference &operator=(const T &v) {
            tree->update(index, v);
            return *this;
        }

        Reference &operator=(const Reference &r) {
            return *this = T(r);
        }
    };

    Reference operator[](Int i) {
        return {this, i};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
    T get_all() const {
        return product(root);
    }
};

template <class T, class Op> auto initDynamicSegmentTree(Int n, Op op, T e) {
    return DynamicSegmentTree<T>(n, op, e);
}

template <class Op, class T> auto newDynamicSegWith(Int n, Op op, T e) {
    return initDynamicSegmentTree(n, op, e);
}

template <class T> Int len(const DynamicSegmentTree<T> &s) {
    return s.len();
}

template <class T> Int node_count(const DynamicSegmentTree<T> &s) {
    return s.node_count();
}

template <class T> void update(DynamicSegmentTree<T> &s, Int i, T value) {
    s.update(i, value);
}

template <class T> T get(const DynamicSegmentTree<T> &s, Int l, Int r) {
    return s.get(l, r);
}

template <class T, class L, class R> T get(const DynamicSegmentTree<T> &s, ClosedSlice<L, R> r) {
    return s.get(r);
}

template <class T> T get_all(const DynamicSegmentTree<T> &s) {
    return s.get_all();
}
}
