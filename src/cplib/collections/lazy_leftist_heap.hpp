#pragma once
#include <cplib/common.hpp>

namespace cplib {
template <class K, class V> class LazyLeftistHeapPool {
    struct Node {
        K key, lazy;
        Int left = -1, right = -1, rank = 1;
        V value;
    };
    std::vector<Node> nodes_;
    K zero_;
    void propagate(Int root) {
        addAll(nodes_[root].left, nodes_[root].lazy);
        addAll(nodes_[root].right, nodes_[root].lazy);
        nodes_[root].lazy = zero_;
    }
    Int rank(Int root) const {
        return root < 0 ? 0 : nodes_[root].rank;
    }
public:
    // 空の根を-1とするヒーププールを作る。zeroは加法単位元。時間・領域 O(capacity)。
    explicit LazyLeftistHeapPool(Int capacity = 0, K zero = K{}) : zero_(zero) {
        assert(capacity >= 0);
        nodes_.reserve(capacity);
    }
    // 一要素の根を作る。同値なら挿入順に取り出す。償却 O(1)。
    Int singleton(K key, V value) {
        Int index = nodes_.size();
        nodes_.push_back({key, zero_, -1, -1, 1, value});
        return index;
    }
    // 全キーにdeltaを加える。空は無視する。計算量 O(1)。
    void addAll(Int root, K delta) {
        if (root >= 0) {
            nodes_[root].key += delta;
            nodes_[root].lazy += delta;
        }
    }
    // 同一プールの互いに素なヒープを破壊的に併合する。時間・再帰領域 O(log N)。
    Int meld(Int a, Int b) {
        if (a < 0)
            return b;
        if (b < 0)
            return a;
        if (nodes_[b].key < nodes_[a].key || (nodes_[a].key == nodes_[b].key && b < a))
            std::swap(a, b);
        propagate(a);
        nodes_[a].right = meld(nodes_[a].right, b);
        if (rank(nodes_[a].left) < rank(nodes_[a].right))
            std::swap(nodes_[a].left, nodes_[a].right);
        nodes_[a].rank = rank(nodes_[a].right) + 1;
        return a;
    }
    // 空でないヒープの最小キーと値を返す。計算量 O(1)。
    std::pair<K, V> top(Int root) const {
        assert(root >= 0);
        return {nodes_[root].key, nodes_[root].value};
    }
    // 最小要素を除いた根を返す。削除前の根は再利用不可。時間 O(log N)。
    Int pop(Int root) {
        assert(root >= 0);
        propagate(root);
        return meld(nodes_[root].left, nodes_[root].right);
    }
};
// ヒーププールを初期化する。時間・領域 O(capacity)。
template <class K, class V>
auto initLazyLeftistHeapPool(Int capacity = 0, K zero = K{}) {
    return LazyLeftistHeapPool<K, V>(capacity, zero);
}
}
