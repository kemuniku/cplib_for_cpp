#pragma once
#include <cplib/common.hpp>
#include <deque>
#include <ranges>
#include <type_traits>

namespace cplib {
template <class T> class MergeSortTree {
    // 子を左、右の順に格納してから親を置く。deque への末尾追加は子への参照を保つ。
    std::deque<T> nodes_;

    template <class Leaves, class Merge>
    const T &build(const Leaves &leaves, Merge &merge, Int first, Int last) {
        if (last - first == 1) {
            return nodes_.emplace_back(leaves[first]);
        }
        const Int mid = std::midpoint(first, last);
        const T &left = build(leaves, merge, first, mid);
        const T &right = build(leaves, merge, mid, last);
        return nodes_.emplace_back(merge(left, right));
    }

    template <class Tree, class Visit>
    static void visitRange(Tree &tree, Int l, Int r, Visit &visit) {
        assert(0 <= l && l <= r && r <= tree.size());
        if (l == r) {
            return;
        }

        auto traverse = [&](auto &self, Int node, Int first, Int last) -> void {
            if (l <= first && last <= r) {
                visit(tree.nodes_[node]);
                return;
            }
            const Int mid = std::midpoint(first, last);
            // 親の直前にある右部分木は、葉が k 個なら全体で 2k - 1 ノードを占める。
            if (l < mid) {
                const Int left = node - 2 * (last - mid);
                self(self, left, first, mid);
            }
            if (mid < r) {
                self(self, node - 1, mid, last);
            }
        };
        traverse(traverse, static_cast<Int>(tree.nodes_.size()) - 1, 0, tree.size());
    }

public:
    // 葉をコピーし、各内部ノードを merge(左の子, 右の子) で構築する。
    // 葉 N 個のコピーと max(N - 1, 0) 回のマージを行う。葉のコピーが O(1) で
    // マージが区間長に比例する場合、構築時間は O(1 + N log(N + 1))。
    template <class Leaves, class Merge>
    MergeSortTree(const Leaves &leaves, Merge &&merge) {
        const Int n = std::size(leaves);
        if (n > 0) {
            build(leaves, merge, 0, n);
        }
    }

    // 葉の個数を O(1) で返す。空でない木のノード数は 2N - 1 個。
    Int size() const {
        return static_cast<Int>((nodes_.size() + 1) / 2);
    }

    // 半開区間 [l, r) を構成するノードの値を左から順に参照渡しする。
    // 値を書き換えても他のノードには反映しない。訪問処理を除いて O(log(N + 1))、空区間は O(1)。
    template <class Visit> void forEachNode(Int l, Int r, Visit &&visit) {
        visitRange(*this, l, r, visit);
    }

    // 半開区間 [l, r) の各ノードを左から順に読み取り専用で訪問する。
    // 訪問処理を除いて O(log(N + 1))、空区間は O(1) で何も呼ばない。
    template <class Visit> void forEachNode(Int l, Int r, Visit &&visit) const {
        visitRange(*this, l, r, visit);
    }
};

// 配列の要素型を葉の型として推論する。実行時の費用はない。
template <class Leaves, class Merge>
MergeSortTree(const Leaves &, Merge &&) -> MergeSortTree<std::ranges::range_value_t<Leaves>>;

// N 個の葉から木を構築する。型を省略すると葉の要素型を使う。
// 葉のコピー N 回と左右の子のマージ max(N - 1, 0) 回を行う。
template <class T = void, class Leaves, class Merge>
auto initMergeSortTree(const Leaves &leaves, Merge &&merge) {
    using Value = std::conditional_t<std::is_void_v<T>, std::ranges::range_value_t<Leaves>, T>;
    return MergeSortTree<Value>(leaves, std::forward<Merge>(merge));
}

// 半開区間 [l, r) の各ノードを左から順に参照渡しする。訪問処理を除いて O(log(N + 1))、空区間は O(1)。
template <class T, class Visit>
void forEachNode(MergeSortTree<T> &tree, Int l, Int r, Visit &&visit) {
    tree.forEachNode(l, r, std::forward<Visit>(visit));
}

// 半開区間 [l, r) の各ノードを左から順に読み取り専用で訪問する。訪問処理を除いて O(log(N + 1))、空区間は O(1)。
template <class T, class Visit>
void forEachNode(const MergeSortTree<T> &tree, Int l, Int r, Visit &&visit) {
    tree.forEachNode(l, r, std::forward<Visit>(visit));
}
}
