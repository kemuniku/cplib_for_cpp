#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <optional>
#include <sstream>
#include <cplib/utils/backwards_index.hpp>

namespace cplib::raw_ptr_avl {
// 生ポインタ版。ノードの所有権は集合のプールが管理する。
template <class K> struct AvlTreeNode {
    using Ptr = AvlTreeNode *;
    Ptr l = nullptr, r = nullptr, p = nullptr;
    Int h = 1, len = 1;
    K key;

    // keyを指定して孤立したノードを構築する。O(1)とkeyのムーブに要する時間。
    explicit AvlTreeNode(K value) : key(std::move(value)) {
    }
};
template <class K> using AvlNodePtr = AvlTreeNode<K> *;

namespace detail {
template <class K> void avl_update(AvlNodePtr<K> n) {
    n->h = 1;
    n->len = 1;
    if (n->l) {
        n->h = std::max(n->h, n->l->h + 1);
        n->len += n->l->len;
    }
    if (n->r) {
        n->h = std::max(n->h, n->r->h + 1);
        n->len += n->r->len;
    }
}

template <class K> void avl_children(AvlNodePtr<K> n, AvlNodePtr<K> l, AvlNodePtr<K> r) {
    n->l = l;
    if (l)
        l->p = n;
    n->r = r;
    if (r)
        r->p = n;
    avl_update(n);
}

template <class K> AvlNodePtr<K> avl_rebalance(AvlNodePtr<K> n) {
    auto l = n->l, r = n->r;
    Int lh = l ? l->h : 0, rh = r ? r->h : 0;
    if (lh + 1 < rh) {
        auto rl = r->l, rr = r->r;
        if ((rl ? rl->h : 0) <= (rr ? rr->h : 0)) {
            r->p = n->p;
            avl_children(n, l, rl);
            avl_children(r, n, rr);
            return r;
        } else {
            rl->p = n->p;
            auto rll = rl->l, rlr = rl->r;
            avl_children(n, l, rll);
            avl_children(r, rlr, rr);
            avl_children(rl, n, r);
            return rl;
        }
    }
    if (rh + 1 < lh) {
        auto ll = l->l, lr = l->r;
        if ((lr ? lr->h : 0) <= (ll ? ll->h : 0)) {
            l->p = n->p;
            avl_children(n, lr, r);
            avl_children(l, ll, n);
            return l;
        } else {
            lr->p = n->p;
            auto lrl = lr->l, lrr = lr->r;
            avl_children(n, lrr, r);
            avl_children(l, ll, lrl);
            avl_children(lr, l, n);
            return lr;
        }
    }
    avl_update(n);
    return n;
}

template <class K> AvlNodePtr<K> avl_to_root(AvlNodePtr<K> n) {
    while (auto p = n->p) {
        if (p->l == n)
            p->l = avl_rebalance(n);
        else
            p->r = avl_rebalance(n);
        n = p;
    }
    return avl_rebalance(n);
}

template <class K> auto avl_search(AvlNodePtr<K> n, const K &key, bool strict) {
    AvlNodePtr<K> l = nullptr, r = nullptr;
    while (n) {
        if (strict ? key < n->key : !(n->key < key)) {
            r = n;
            n = n->l;
        } else {
            l = n;
            n = n->r;
        }
    }
    return std::pair{l, r};
}
}

// ノードnから親をたどって根を返す。木の高さをhとしてO(h)。
template <class K> AvlNodePtr<K> rootOf(AvlNodePtr<K> n) {
    while (auto p = n->p)
        n = p;
    return n;
}

// key未満の最大ノードとkey以上の最小ノードの組を返す。木の高さをhとしてO(h)。
template <class K> auto lower_bound_node(AvlNodePtr<K> n, const K &key) {
    return detail::avl_search(n, key, false);
}

// key以下の最大ノードとkeyより大きい最小ノードの組を返す。木の高さをhとしてO(h)。
template <class K> auto upper_bound_node(AvlNodePtr<K> n, const K &key) {
    return detail::avl_search(n, key, true);
}

// 探索済みの前後ノードを使ってxを挿入し、新しい根を返す。要素数をnとしてO(log(n+1))。
template <class K>
AvlNodePtr<K> insertAfterLowerBound(AvlNodePtr<K> n, AvlNodePtr<K> x, AvlNodePtr<K> l,
                                    AvlNodePtr<K> r) {
    if (!n)
        return x;
    if (l && !l->r) {
        detail::avl_children(l, l->l, x);
        return detail::avl_to_root(l);
    }
    detail::avl_children(r, x, r->r);
    return detail::avl_to_root(r);
}

// 根nのAVL木にノードxを挿入し、新しい根を返す。要素数をnとしてO(log(n+1))。
template <class K> AvlNodePtr<K> insert(AvlNodePtr<K> n, AvlNodePtr<K> x) {
    auto [l, r] = lower_bound_node(n, x->key);
    return insertAfterLowerBound(n, x, l, r);
}

// ノードxを削除し、新しい根を返す。nxtにはxの次ノードを渡す。要素数をnとしてO(log(n+1))。
template <class K> AvlNodePtr<K> erase(AvlNodePtr<K>, AvlNodePtr<K> x, AvlNodePtr<K> nxt) {
    auto xp = x->p;
    AvlNodePtr<K> result = nullptr;
    if (!x->r) {
        auto xl = x->l;
        if (xl)
            xl->p = xp;
        if (xp) {
            if (xp->l == x)
                xp->l = xl;
            else
                xp->r = xl;
        }
        result = xp ? detail::avl_to_root(xp) : xl;
    } else {
        auto nxtp = nxt->p, nxtr = nxt->r;
        if (xp) {
            if (xp->l == x)
                xp->l = nxt;
            else
                xp->r = nxt;
        }
        nxt->p = xp;
        nxt->l = x->l;
        if (nxt->l)
            nxt->l->p = nxt;
        if (x->r == nxt) {
            detail::avl_update(nxt);
            result = detail::avl_to_root(nxt);
        } else {
            if (nxtp->l == nxt)
                nxtp->l = nxtr;
            else
                nxtp->r = nxtr;
            if (nxtr)
                nxtr->p = nxtp;
            nxt->r = x->r;
            nxt->r->p = nxt;
            detail::avl_update(nxt);
            result = detail::avl_to_root(nxtp);
        }
    }
    x->l = x->r = x->p = nullptr;
    detail::avl_update(x);
    return result;
}

// 昇順で次のノードを返す。木の高さをhとしてO(h)。
template <class K> AvlNodePtr<K> next(AvlNodePtr<K> n) {
    if (n->r) {
        n = n->r;
        while (n->l)
            n = n->l;
        return n;
    }
    while (auto p = n->p) {
        if (p->r != n)
            return p;
        n = p;
    }
    return {};
}

// 昇順で前のノードを返す。木の高さをhとしてO(h)。
template <class K> AvlNodePtr<K> prev(AvlNodePtr<K> n) {
    if (n->l) {
        n = n->l;
        while (n->r)
            n = n->r;
        return n;
    }
    while (auto p = n->p) {
        if (p->l != n)
            return p;
        n = p;
    }
    return {};
}

// 部分木nの昇順でidx番目のノードを返す。木の高さをhとしてO(h)。
template <class K> AvlNodePtr<K> get(AvlNodePtr<K> n, Int idx) {
    assert(idx >= 0);
    if (!n || idx >= n->len)
        return {};
    for (;;) {
        Int left = n->l ? n->l->len : 0;
        if (left == idx)
            return n;
        if (left < idx) {
            idx -= left + 1;
            n = n->r;
        } else
            n = n->l;
    }
}

// ノードnの木全体での昇順順位を返す。木の高さをhとしてO(h)。
template <class K> Int index(AvlNodePtr<K> n) {
    if (!n)
        return 0;
    Int result = n->l ? n->l->len : 0;
    while (auto p = n->p) {
        if (p->r == n)
            result += 1 + (p->l ? p->l->len : 0);
        n = p;
    }
    return result;
}

// 順位付き AVL 集合。挿入・削除・順位取得 O(log N)、列挙 O(N)。
template <class T, bool Multi> class BasicAvlSortedSet {
    std::vector<std::unique_ptr<AvlTreeNode<T>>> nodes_;

    AvlNodePtr<T> newNode(const T &x) {
        nodes_.push_back(std::make_unique<AvlTreeNode<T>>(x));
        return nodes_.back().get();
    }

    AvlNodePtr<T> clone(AvlNodePtr<T> n, AvlNodePtr<T> parent = nullptr) {
        if (!n)
            return nullptr;
        auto out = newNode(n->key);
        out->h = n->h;
        out->len = n->len;
        out->p = parent;
        out->l = clone(n->l, out);
        out->r = clone(n->r, out);
        return out;
    }

    AvlNodePtr<T> buildBalanced(const std::vector<T> &a, Int l, Int r,
                                AvlNodePtr<T> parent = nullptr) {
        if (l >= r)
            return nullptr;
        Int m = (l + r) / 2;
        auto out = newNode(a[m]);
        out->p = parent;
        out->l = buildBalanced(a, l, m, out);
        out->r = buildBalanced(a, m + 1, r, out);
        detail::avl_update(out);
        return out;
    }

    static std::optional<T> value(AvlNodePtr<T> n) {
        return n ? std::optional<T>(n->key) : std::nullopt;
    }

public:
    AvlNodePtr<T> root = nullptr;
    // 空の集合を構築する。O(1)。
    BasicAvlSortedSet() = default;

    // 配列の要素を挿入して構築する。長さをnとしてO(n log(n+1))。
    explicit BasicAvlSortedSet(std::vector<T> v) {
        std::sort(v.begin(), v.end());
        if constexpr (!Multi)
            v.erase(std::unique(v.begin(), v.end()), v.end());
        root = buildBalanced(v, 0, v.size());
    }

    // 全要素を複製する。要素数をnとしてO(n)。
    BasicAvlSortedSet(const BasicAvlSortedSet &other) : root(clone(other.root)) {
    }

    // 集合を複製して代入する。代入元の要素数をn、代入先の確保済みノード数をMとしてO(n+M)。
    BasicAvlSortedSet &operator=(const BasicAvlSortedSet &other) {
        if (this != &other) {
            nodes_.clear();
            root = clone(other.root);
        }
        return *this;
    }

    // 集合をムーブ構築する。O(1)。
    BasicAvlSortedSet(BasicAvlSortedSet &&other) noexcept
        : nodes_(std::move(other.nodes_)), root(std::exchange(other.root, nullptr)) {
    }

    // 集合をムーブ代入する。代入先の確保済みノード数をMとしてO(M)。
    BasicAvlSortedSet &operator=(BasicAvlSortedSet &&other) noexcept {
        if (this != &other) {
            nodes_ = std::move(other.nodes_);
            root = std::exchange(other.root, nullptr);
        }
        return *this;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return root ? root->len : 0;
    }

    // x未満の要素数を返す。O(log(n+1))。
    Int lowerBound(const T &x) const {
        auto n = lower_bound_node(root, x).second;
        return n ? cplib::raw_ptr_avl::index(n) : size();
    }

    // x以下の要素数を返す。O(log(n+1))。
    Int upperBound(const T &x) const {
        auto n = upper_bound_node(root, x).second;
        return n ? cplib::raw_ptr_avl::index(n) : size();
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
        return n && n->key == x;
    }

    // xを追加する。重複を許さない集合で既存ならfalse。O(log(n+1))。
    bool incl(const T &x) {
        auto [l, r] = lower_bound_node(root, x);
        if constexpr (!Multi) {
            if (r && r->key == x)
                return false;
        }
        root = insertAfterLowerBound(root, newNode(x), l, r);
        return true;
    }

    // xを一つ削除し、削除できたか返す。O(log(n+1))。
    bool excl(const T &x) {
        auto n = lower_bound_node(root, x).second;
        if (!n || !(n->key == x))
            return false;
        root = cplib::raw_ptr_avl::erase(root, n, cplib::raw_ptr_avl::next(n));
        return true;
    }

    // 昇順でi番目の要素を返す。O(log(n+1))。
    T operator[](Int i) const {
        assert(i >= 0 && i < size());
        return get(root, i)->key;
    }

    // 昇順の末尾からi番目の要素を返す。O(log(n+1))。
    T operator[](BackwardsIndex i) const {
        return (*this)[size() - i.value];
    }

    // 昇順でi番目の要素を削除して返す。負のiは末尾から数える。O(log(n+1))。
    T pop(Int i = -1) {
        if (i < 0)
            i += size();
        assert(i >= 0 && i < size());
        auto n = get(root, i);
        T out = n->key;
        root = cplib::raw_ptr_avl::erase(root, n, cplib::raw_ptr_avl::next(n));
        return out;
    }

    struct Iterator {
        AvlNodePtr<T> node = nullptr;

        // 指している要素への参照を返す。O(1)。
        const T &operator*() const {
            return node->key;
        }

        // 次の要素へ進む。最悪O(log(n+1))、全列挙では合計O(n)。
        Iterator &operator++() {
            node = cplib::raw_ptr_avl::next(node);
            return *this;
        }

        // 同じ位置を指すか比較する。O(1)。
        bool operator==(const Iterator &) const = default;
    };

    // 先頭を指すイテレータを返す。O(log(n+1))。
    Iterator begin() const {
        auto n = root;
        if (n)
            while (n->l)
                n = n->l;
        return {n};
    }

    // 末尾の次を指すイテレータを返す。O(1)。
    Iterator end() const {
        return {};
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
