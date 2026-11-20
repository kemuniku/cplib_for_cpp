#pragma once
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

namespace cplib::detail {

// 走査線上の辺の順序を保持する永続 AVL 木。比較関数は更新ごとに指定する。
class PersistentEdgeOrder {
public:
    struct Node {
        int edge, left, right, height;
    };

private:
    std::vector<Node> nodes_;
    std::pair<int, int> remove_min(int root) {
        Node a = nodes_[root];
        if (a.left < 0)
            return {a.edge, a.right};
        auto [edge, left] = remove_min(a.left);
        return {edge, balance(make(a.edge, left, a.right))};
    }
    int height(int node) const { return node < 0 ? 0 : nodes_[node].height; }
    int make(int edge, int left, int right) {
        nodes_.push_back({edge, left, right, 1 + std::max(height(left), height(right))});
        return int(nodes_.size()) - 1;
    }
    int rotate_left(int root) {
        Node a = nodes_[root], b = nodes_[a.right];
        int child = make(a.edge, a.left, b.left);
        return make(b.edge, child, b.right);
    }
    int rotate_right(int root) {
        Node a = nodes_[root], b = nodes_[a.left];
        int child = make(a.edge, b.right, a.right);
        return make(b.edge, b.left, child);
    }
    int balance(int root) {
        Node a = nodes_[root];
        if (height(a.left) - height(a.right) > 1) {
            Node b = nodes_[a.left];
            if (height(b.left) < height(b.right))
                root = make(a.edge, rotate_left(a.left), a.right);
            return rotate_right(root);
        }
        if (height(a.right) - height(a.left) > 1) {
            Node b = nodes_[a.right];
            if (height(b.right) < height(b.left))
                root = make(a.edge, a.left, rotate_right(a.right));
            return rotate_left(root);
        }
        return root;
    }

public:

    const Node &operator[](int node) const { return nodes_[node]; }

    template <class Compare> int insert(int root, int edge, const Compare &compare) {
        if (root < 0)
            return make(edge, -1, -1);
        Node a = nodes_[root];
        int cmp = compare(edge, a.edge);
        assert(cmp != 0);
        if (cmp < 0)
            a.left = insert(a.left, edge, compare);
        else
            a.right = insert(a.right, edge, compare);
        return balance(make(a.edge, a.left, a.right));
    }
    template <class Compare> int erase(int root, int edge, const Compare &compare) {
        assert(root >= 0);
        Node a = nodes_[root];
        int cmp = compare(edge, a.edge);
        if (cmp < 0)
            return balance(make(a.edge, erase(a.left, edge, compare), a.right));
        if (cmp > 0)
            return balance(make(a.edge, a.left, erase(a.right, edge, compare)));
        if (a.left < 0)
            return a.right;
        if (a.right < 0)
            return a.left;
        auto [next, right] = remove_min(a.right);
        return balance(make(next, a.left, right));
    }
};
}
