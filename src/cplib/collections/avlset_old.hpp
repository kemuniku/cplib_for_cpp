#pragma once
#include <cplib/collections/avltreenode_old.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <optional>
#include <sstream>

namespace cplib {
// 順位付き AVL 集合。挿入・削除・順位取得 O(log N)、列挙 O(N)。
template <class T, bool Multi> class BasicAvlSortedSet {
    static AvlNodePtr<T> clone(AvlNodePtr<T> n, AvlNodePtr<T> parent = get_avltree_nilnode<T>()) {
        if (n == get_avltree_nilnode<T>())
            return get_avltree_nilnode<T>();
        auto out = std::make_shared<AvlTreeNode<T>>(n->key);
        out->h = n->h;
        out->len = n->len;
        out->p = parent;
        out->l = clone(n->l, out);
        out->r = clone(n->r, out);
        return out;
    }

    static std::optional<T> value(AvlNodePtr<T> n) {
        return n != get_avltree_nilnode<T>() ? std::optional<T>(n->key) : std::nullopt;
    }

public:
    AvlNodePtr<T> root = get_avltree_nilnode<T>();
    // 空の集合を構築する。O(1)。
    BasicAvlSortedSet() = default;

    // 配列の要素を挿入して構築する。長さをnとしてO(n log(n+1))。
    explicit BasicAvlSortedSet(const std::vector<T> &v) {
        for (const T &x : v)
            incl(x);
    }

    // 全要素を複製する。要素数をnとしてO(n)。
    BasicAvlSortedSet(const BasicAvlSortedSet &other) : root(clone(other.root)) {
    }

    // 集合を複製して代入する。両者の要素数の合計をnとしてO(n)。
    BasicAvlSortedSet &operator=(const BasicAvlSortedSet &other) {
        if (this != &other)
            root = clone(other.root);
        return *this;
    }

    // 集合をムーブ構築する。O(1)。
    BasicAvlSortedSet(BasicAvlSortedSet &&) = default;
    // 集合をムーブ代入する。代入先の要素数をnとしてO(n)。
    BasicAvlSortedSet &operator=(BasicAvlSortedSet &&) = default;

    // 要素数を返す。O(1)。
    Int len() const {
        return root->len;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return len();
    }

    // x未満の要素数を返す。O(log(n+1))。
    Int lowerBound(const T &x) const {
        auto n = lower_bound_node(root, x).second;
        return n != get_avltree_nilnode<T>() ? cplib::index(n) : len();
    }

    // x以下の要素数を返す。O(log(n+1))。
    Int upperBound(const T &x) const {
        auto n = upper_bound_node(root, x).second;
        return n != get_avltree_nilnode<T>() ? cplib::index(n) : len();
    }

    // x未満の要素数を返す。O(log(n+1))。
    Int index(const T &x) const {
        return lowerBound(x);
    }

    // x以下の要素数を返す。O(log(n+1))。
    Int index_right(const T &x) const {
        return upperBound(x);
    }

    // xの個数を返す。O(log(n+1))。
    Int count(const T &x) const {
        return upperBound(x) - lowerBound(x);
    }

    // x未満の最大値を返す。存在しなければ空。O(log(n+1))。
    auto lt(const T &x) const {
        return value(lower_bound_node(root, x).first);
    }

    // x以下の最大値を返す。存在しなければ空。O(log(n+1))。
    auto le(const T &x) const {
        return value(upper_bound_node(root, x).first);
    }

    // xより大きい最小値を返す。存在しなければ空。O(log(n+1))。
    auto gt(const T &x) const {
        return value(upper_bound_node(root, x).second);
    }

    // x以上の最小値を返す。存在しなければ空。O(log(n+1))。
    auto ge(const T &x) const {
        return value(lower_bound_node(root, x).second);
    }

    // xが存在するか返す。O(log(n+1))。
    bool contains(const T &x) const {
        auto n = lower_bound_node(root, x).second;
        return n != get_avltree_nilnode<T>() && n->key == x;
    }

    // xを追加する。重複を許さない集合で既存ならfalse。O(log(n+1))。
    bool incl(const T &x) {
        if constexpr (!Multi) {
            if (contains(x))
                return false;
        }
        root = cplib::insert(root, std::make_shared<AvlTreeNode<T>>(x));
        return true;
    }

    // xを一つ削除し、削除できたか返す。O(log(n+1))。
    bool excl(const T &x) {
        auto n = lower_bound_node(root, x).second;
        if (n == get_avltree_nilnode<T>() || !(n->key == x))
            return false;
        root = cplib::erase(root, n, cplib::next(n));
        return true;
    }

    // 昇順でi番目の要素を返す。O(log(n+1))。
    T operator[](Int i) const {
        assert(i >= 0 && i < len());
        return get(root, i)->key;
    }

    // 昇順の末尾からi番目の要素を返す。O(log(n+1))。
    T operator[](BackwardsIndex i) const {
        return (*this)[len() - i.value];
    }

    // 昇順でi番目の要素を削除して返す。負のiは末尾から数える。O(log(n+1))。
    T pop(Int i = -1) {
        if (i < 0)
            i += len();
        assert(i >= 0 && i < len());
        auto n = get(root, i);
        T out = n->key;
        root = cplib::erase(root, n, cplib::next(n));
        return out;
    }

    struct Iterator {
        AvlNodePtr<T> node;

        // 指している要素への参照を返す。O(1)。
        const T &operator*() const {
            return node->key;
        }

        // 次の要素へ進む。最悪O(log(n+1))、全列挙では合計O(n)。
        Iterator &operator++() {
            node = cplib::next(node);
            return *this;
        }

        // 同じ位置を指すか比較する。O(1)。
        bool operator==(const Iterator &) const = default;
    };

    // 先頭を指すイテレータを返す。O(log(n+1))。
    Iterator begin() const {
        auto n = root;
        if (n != get_avltree_nilnode<T>())
            while (n->l != get_avltree_nilnode<T>())
                n = n->l;
        return {n};
    }

    // 末尾の次を指すイテレータを返す。O(1)。
    Iterator end() const {
        return {get_avltree_nilnode<T>()};
    }

    // 列挙用に自身への参照を返す。O(1)。
    const BasicAvlSortedSet &items() const {
        return *this;
    }

    // 要素を昇順に空白区切りで文字列化する。O(n+出力長)と各要素の文字列化の計算量。
    std::string to_string() const {
        std::ostringstream out;
        bool first = true;
        for (const T &x : *this) {
            if (!first)
                out << ' ';
            first = false;
            out << x;
        }
        return out.str();
    }
};
template <class T> using AvlSortedMultiSet = BasicAvlSortedSet<T, true>;
template <class T> using AVLSortedMultiSet = AvlSortedMultiSet<T>;
template <class T> using AvlSortedSet = BasicAvlSortedSet<T, false>;
template <class T> using AVLSortedSet = AvlSortedSet<T>;

// 配列から重複を許す集合を構築する。長さをnとしてO(n log(n+1))。
template <class T> auto initAvlSortedMultiSet(const std::vector<T> &v = {}) {
    return AvlSortedMultiSet<T>(v);
}

// 配列から重複を許さない集合を構築する。長さをnとしてO(n log(n+1))。
template <class T> auto initAvlSortedSet(const std::vector<T> &v = {}) {
    return AvlSortedSet<T>(v);
}
}
