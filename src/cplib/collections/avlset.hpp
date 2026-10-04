#pragma once
#include <cplib/collections/avltreenode.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <optional>
#include <sstream>

namespace cplib {
// 順位付き AVL 集合。挿入・削除・順位取得 O(log N)、列挙 O(N)。
template <class T, bool Multi> class BasicAvlSortedSet {
    static AvlNodePtr<T> clone(AvlNodePtr<T> n, AvlNodePtr<T> parent = {}) {
        if (!n)
            return {};
        auto out = std::make_shared<AvlTreeNode<T>>(n->key);
        out->h = n->h;
        out->len = n->len;
        out->p = parent;
        out->l = clone(n->l, out);
        out->r = clone(n->r, out);
        return out;
    }

    static std::optional<T> value(AvlNodePtr<T> n) {
        return n ? std::optional<T>(n->key) : std::nullopt;
    }

public:
    AvlNodePtr<T> root;
    BasicAvlSortedSet() = default;

    explicit BasicAvlSortedSet(const std::vector<T> &v) {
        for (const T &x : v)
            incl(x);
    }

    BasicAvlSortedSet(const BasicAvlSortedSet &other) : root(clone(other.root)) {
    }

    BasicAvlSortedSet &operator=(const BasicAvlSortedSet &other) {
        if (this != &other)
            root = clone(other.root);
        return *this;
    }

    BasicAvlSortedSet(BasicAvlSortedSet &&) = default;
    BasicAvlSortedSet &operator=(BasicAvlSortedSet &&) = default;

    Int len() const {
        return root ? root->len : 0;
    }

    Int size() const {
        return len();
    }

    Int lowerBound(const T &x) const {
        auto n = lower_bound_node(root, x).second;
        return n ? cplib::index(n) : len();
    }

    Int upperBound(const T &x) const {
        auto n = upper_bound_node(root, x).second;
        return n ? cplib::index(n) : len();
    }

    Int index(const T &x) const {
        return lowerBound(x);
    }

    Int index_right(const T &x) const {
        return upperBound(x);
    }

    Int count(const T &x) const {
        return upperBound(x) - lowerBound(x);
    }

    auto lt(const T &x) const {
        return value(lower_bound_node(root, x).first);
    }

    auto le(const T &x) const {
        return value(upper_bound_node(root, x).first);
    }

    auto gt(const T &x) const {
        return value(upper_bound_node(root, x).second);
    }

    auto ge(const T &x) const {
        return value(lower_bound_node(root, x).second);
    }

    bool contains(const T &x) const {
        auto n = lower_bound_node(root, x).second;
        return n && n->key == x;
    }

    bool incl(const T &x) {
        if constexpr (!Multi) {
            if (contains(x))
                return false;
        }
        root = cplib::insert(root, std::make_shared<AvlTreeNode<T>>(x));
        return true;
    }

    bool excl(const T &x) {
        auto n = lower_bound_node(root, x).second;
        if (!n || !(n->key == x))
            return false;
        root = cplib::erase(root, n, cplib::next(n));
        return true;
    }

    T operator[](Int i) const {
        assert(i >= 0 && i < len());
        return get(root, i)->key;
    }

    T operator[](BackwardsIndex i) const {
        return (*this)[len() - i.value];
    }

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

        const T &operator*() const {
            return node->key;
        }

        Iterator &operator++() {
            node = cplib::next(node);
            return *this;
        }

        bool operator==(const Iterator &) const = default;
    };

    Iterator begin() const {
        auto n = root;
        if (n)
            while (n->l)
                n = n->l;
        return {n};
    }

    Iterator end() const {
        return {};
    }

    const BasicAvlSortedSet &items() const {
        return *this;
    }

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

template <class T> auto initAvlSortedMultiSet(const std::vector<T> &v = {}) {
    return AvlSortedMultiSet<T>(v);
}

template <class T> auto initAvlSortedSet(const std::vector<T> &v = {}) {
    return AvlSortedSet<T>(v);
}
}
