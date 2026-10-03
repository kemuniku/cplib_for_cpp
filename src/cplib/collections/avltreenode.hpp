#pragma once
#include <cplib/common.hpp>
#include <memory>

namespace cplib {
// 子が所有し、親は弱参照。回転・削除後も外部が保持したノード参照は有効。
template <class K> struct AvlTreeNode {
    using Ptr = std::shared_ptr<AvlTreeNode>;
    Ptr l, r;
    std::weak_ptr<AvlTreeNode> p;
    Int h = 1, len = 1;
    K key;

    explicit AvlTreeNode(K value) : key(std::move(value)) {
    }
};
template <class K> using AvlNodePtr = std::shared_ptr<AvlTreeNode<K>>;

namespace detail {
template <class K> void avl_update(AvlNodePtr<K> n) {
    n->h = 0;
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
    while (auto p = n->p.lock()) {
        if (p->l == n)
            p->l = avl_rebalance(n);
        else
            p->r = avl_rebalance(n);
        n = p;
    }
    return avl_rebalance(n);
}

template <class K> auto avl_search(AvlNodePtr<K> n, const K &key, bool strict) {
    AvlNodePtr<K> l, r;
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

template <class K> AvlNodePtr<K> rootOf(AvlNodePtr<K> n) {
    while (auto p = n->p.lock())
        n = p;
    return n;
}

template <class K> auto lower_bound_node(AvlNodePtr<K> n, const K &key) {
    return detail::avl_search(n, key, false);
}

template <class K> auto upper_bound_node(AvlNodePtr<K> n, const K &key) {
    return detail::avl_search(n, key, true);
}

template <class K> AvlNodePtr<K> insert(AvlNodePtr<K> n, AvlNodePtr<K> x) {
    if (!n)
        return x;
    auto [l, r] = lower_bound_node(n, x->key);
    if (l && !l->r) {
        detail::avl_children(l, l->l, x);
        return detail::avl_to_root(l);
    }
    detail::avl_children(r, x, r->r);
    return detail::avl_to_root(r);
}

template <class K> AvlNodePtr<K> erase(AvlNodePtr<K>, AvlNodePtr<K> x, AvlNodePtr<K> nxt) {
    auto xp = x->p.lock();
    AvlNodePtr<K> result;
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
        auto nxtp = nxt->p.lock(), nxtr = nxt->r;
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
    x->l.reset();
    x->r.reset();
    x->p.reset();
    detail::avl_update(x);
    return result;
}

template <class K> AvlNodePtr<K> next(AvlNodePtr<K> n) {
    if (n->r) {
        n = n->r;
        while (n->l)
            n = n->l;
        return n;
    }
    while (auto p = n->p.lock()) {
        if (p->r != n)
            return p;
        n = p;
    }
    return {};
}

template <class K> AvlNodePtr<K> prev(AvlNodePtr<K> n) {
    if (n->l) {
        n = n->l;
        while (n->r)
            n = n->r;
        return n;
    }
    while (auto p = n->p.lock()) {
        if (p->l != n)
            return p;
        n = p;
    }
    return {};
}

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

template <class K> Int index(AvlNodePtr<K> n) {
    if (!n)
        return 0;
    Int result = n->l ? n->l->len : 0;
    while (auto p = n->p.lock()) {
        if (p->r == n)
            result += 1 + (p->l ? p->l->len : 0);
        n = p;
    }
    return result;
}
}
