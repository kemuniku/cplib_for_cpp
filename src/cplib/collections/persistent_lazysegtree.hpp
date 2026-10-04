#pragma once
#include <cplib/collections/private/persistent_arena.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <functional>
#include <sstream>

namespace cplib {
// Path copying; operations allocate O(log N) nodes, queries allocate none.
// Independent arenas are merged by size on copy_range.
template <class S, class F> class PersistentLazySegmentTree {
    struct Node {
        S value;
        F lazy;
        bool pending = false;
        Node *left = nullptr;
        Node *right = nullptr;
    };

    detail::PersistentArena<Node> arena;
    Node *root = nullptr;
    Int size_ = 1, n = 0;
    std::function<S(S, S)> merge;
    S e;
    std::function<S(F, S)> mapping;
    std::function<F(F, F)> composition;
    F id;

    Node *make(S v, Node *l = nullptr, Node *r = nullptr) const {
        return arena.make({std::move(v), id, false, l, r});
    }

    Node *applied(Node *node, F f) const {
        auto out = make(mapping(f, node->value), node->left, node->right);
        if (node->left) {
            out->lazy = node->pending ? composition(f, node->lazy) : f;
            out->pending = true;
        }
        return out;
    }

    Node *carried(Node *node, const F &carry, bool pending) const {
        return pending ? applied(node, carry) : node;
    }

    std::pair<F, bool> descend(Node *node, const F &carry, bool pending) const {
        if (!node->pending)
            return {carry, pending};
        return {pending ? composition(carry, node->lazy) : node->lazy, true};
    }

    PersistentLazySegmentTree withRoot(Node *p) const {
        auto out = *this;
        out.root = p;
        return out;
    }

    void check(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= n);
    }

public:
    // vから永続遅延セグ木をO(N)で構築します。comp(f,g)はgの後にfを作用させます。
    // opは結合的でidentityは単位元、mapはopと両立し、identity_actionは恒等作用とします。
    // 各演算はO(1)とし、引数や共有する参照先を変更しないでください。
    // 共有・統合したノード領域は、それを所有するすべての版の破棄時に解放します。
    // 個別の古い版を破棄しても、その領域のノードは回収しません。
    template <class Op, class Map, class Comp>
    PersistentLazySegmentTree(std::span<const S> v, Op op, S identity, Map map, Comp comp,
                              F identity_action)
        : n(v.size()), merge(op), e(identity), mapping(map), composition(comp),
          id(identity_action) {
        while (size_ < n)
            size_ *= 2;
        auto build = [&](auto &&self, Int l, Int r) -> Node * {
            if (r - l == 1)
                return make(l < n ? v[l] : e);
            Int m = (l + r) / 2;
            auto a = self(self, l, m), b = self(self, m, r);
            return make(merge(a->value, b->value), a, b);
        };
        root = build(build, 0, size_);
    }

    template <class Op, class Map, class Comp>
    PersistentLazySegmentTree(Int count, Op op, S identity, Map map, Comp comp, F identity_action)
        : PersistentLazySegmentTree(std::vector<S>(checked_count(count), identity), op, identity,
                                    map, comp, identity_action) {
    }

    static std::size_t checked_count(Int count) {
        assert(count >= 0);
        return count;
    }

    Int len() const {
        return n;
    }

    Int size() const {
        return len();
    }

    // pをvalueに置き換えた新しい版を、時間・追加領域O(log N)で返します。
    PersistentLazySegmentTree update(Int p, S value) const {
        assert(0 <= p && p < n);
        auto dfs = [&](auto &&self, Node *node, Int l, Int r, F carry, bool pending) -> Node * {
            if (r - l == 1)
                return make(value);
            Int m = (l + r) / 2;
            auto [f, has] = descend(node, carry, pending);
            Node *a, *b;
            if (p < m) {
                a = self(self, node->left, l, m, f, has);
                b = carried(node->right, f, has);
            } else {
                a = carried(node->left, f, has);
                b = self(self, node->right, m, r, f, has);
            }
            return make(merge(a->value, b->value), a, b);
        };
        return withRoot(dfs(dfs, root, 0, size_, id, false));
    }

    // 半開区間[ql,qr)にfを作用させた新しい版を、時間・追加領域O(log N)で返します。
    PersistentLazySegmentTree apply(Int ql, Int qr, F f) const {
        check(ql, qr);
        if (ql == qr)
            return *this;
        auto dfs = [&](auto &&self, Node *node, Int l, Int r, F carry, bool pending) -> Node * {
            if (qr <= l || r <= ql)
                return carried(node, carry, pending);
            if (ql <= l && r <= qr)
                return applied(node, pending ? composition(f, carry) : f);
            Int m = (l + r) / 2;
            auto [next, has] = descend(node, carry, pending);
            auto a = self(self, node->left, l, m, next, has),
                 b = self(self, node->right, m, r, next, has);
            return make(merge(a->value, b->value), a, b);
        };
        return withRoot(dfs(dfs, root, 0, size_, id, false));
    }

    // 半開区間[ql,qr)をsourceの同じ区間で置き換えた新しい版を返します。木の操作は時間・追加領域O(log N)。
    // 両方の木は同じ長さ・演算・単位元で構築してください。元の版は変更しません。
    // 独立に構築した木の部分コピーでは領域の所有権も統合します。
    // K個の独立した木の統合に伴う領域参照の移動は全操作でO(K log K)。
    PersistentLazySegmentTree copy_range(const PersistentLazySegmentTree &source, Int ql,
                                         Int qr) const {
        assert(n == source.n);
        check(ql, qr);
        if (ql == qr || root == source.root)
            return *this;
        if (ql == 0 && qr == n) {
            auto out = withRoot(source.root);
            out.arena = source.arena;
            return out;
        }
        auto dfs = [&](auto &&self, Node *dst, Node *src, Int l, Int r, F dc, F sc, bool dp,
                       bool sp) -> Node * {
            if (qr <= l || r <= ql)
                return carried(dst, dc, dp);
            if (ql <= l && r <= qr)
                return source.carried(src, sc, sp);
            Int m = (l + r) / 2;
            auto [df, dh] = descend(dst, dc, dp);
            auto [sf, sh] = source.descend(src, sc, sp);
            auto a = self(self, dst->left, src->left, l, m, df, sf, dh, sh),
                 b = self(self, dst->right, src->right, m, r, df, sf, dh, sh);
            return make(merge(a->value, b->value), a, b);
        };
        auto out = withRoot(dfs(dfs, root, source.root, 0, size_, id, source.id, false, false));
        arena.share(source.arena);
        return out;
    }

    // 半開区間[ql,qr)の積をO(log N)で返します。ノードは変更・複製しません。
    S get(Int ql, Int qr) const {
        check(ql, qr);
        if (ql == qr)
            return e;
        auto dfs = [&](auto &&self, Node *node, Int l, Int r) -> S {
            if (ql <= l && r <= qr)
                return node->value;
            Int m = (l + r) / 2;
            S value = qr <= m ? self(self, node->left, l, m)
                      : m <= ql
                          ? self(self, node->right, m, r)
                          : merge(self(self, node->left, l, m), self(self, node->right, m, r));
            return node->pending ? mapping(node->lazy, value) : value;
        };
        return dfs(dfs, root, 0, size_);
    }

    S query(Int l, Int r) const {
        return get(l, r);
    }

    S get_all() const {
        return root->value;
    }

    S operator[](Int p) const {
        assert(0 <= p && p < n);
        return get(p, p + 1);
    }

    // Assignment is explicit to keep ordinary reads value typed.
    void set(Int p, S v) {
        *this = update(p, std::move(v));
    }

    void set(BackwardsIndex p, S v) {
        set(n - p.value, std::move(v));
    }

    S operator[](BackwardsIndex p) const {
        return (*this)[n - p.value];
    }

    template <class L, class R> S get(ClosedSlice<L, R> s) const {
        return get(resolve_index(n, s.a), resolve_index(n, s.b) + 1);
    }

    template <class L, class R> S operator[](ClosedSlice<L, R> s) const {
        return get(s);
    }

    template <class L, class R> auto apply(ClosedSlice<L, R> s, F f) const {
        return apply(resolve_index(n, s.a), resolve_index(n, s.b) + 1, f);
    }

    template <class L, class R>
    auto copy_range(const PersistentLazySegmentTree &src, ClosedSlice<L, R> s) const {
        return copy_range(src, resolve_index(n, s.a), resolve_index(n, s.b) + 1);
    }

    // pred(get(l,r))を満たす最大のrをO(log N)で返します。ノードは変更・複製しません。
    // predは区間の拡大に対して単調で、単位元に対してtrueを返す必要があります。
    template <class Pred> Int max_right(Int l, Pred pred) const {
        assert(0 <= l && l <= n && pred(e));
        S sm = e;
        auto dfs = [&](auto &&self, Node *node, Int nl, Int nr, F carry) -> Int {
            if (nr <= l || n <= nl)
                return n;
            if (l <= nl && nr <= n) {
                S value = merge(sm, mapping(carry, node->value));
                if (pred(value)) {
                    sm = value;
                    return n;
                }
                if (nr - nl == 1)
                    return nl;
            }
            Int m = (nl + nr) / 2;
            F next = composition(carry, node->lazy);
            Int b = self(self, node->left, nl, m, next);
            return b != n ? b : self(self, node->right, m, nr, next);
        };
        return dfs(dfs, root, 0, size_, id);
    }

    // pred(get(l,r))を満たす最小のlをO(log N)で返します。ノードは変更・複製しません。
    // predは区間の拡大に対して単調で、単位元に対してtrueを返す必要があります。
    template <class Pred> Int min_left(Int r, Pred pred) const {
        assert(0 <= r && r <= n && pred(e));
        S sm = e;
        auto dfs = [&](auto &&self, Node *node, Int nl, Int nr, F carry) -> Int {
            if (r <= nl)
                return 0;
            if (nr <= r) {
                S value = merge(mapping(carry, node->value), sm);
                if (pred(value)) {
                    sm = value;
                    return 0;
                }
                if (nr - nl == 1)
                    return nr;
            }
            Int m = (nl + nr) / 2;
            F next = composition(carry, node->lazy);
            Int b = self(self, node->right, m, nr, next);
            return b ? b : self(self, node->left, nl, m, next);
        };
        return dfs(dfs, root, 0, size_, id);
    }

    std::string str() const {
        std::ostringstream out;
        bool first = true;
        auto visit = [&](auto &&self, Node *node, Int l, Int r, F carry) -> void {
            if (l >= n)
                return;
            if (r - l == 1) {
                if (!first)
                    out << ' ';
                first = false;
                out << mapping(carry, node->value);
                return;
            }
            Int m = (l + r) / 2;
            F next = composition(carry, node->lazy);
            self(self, node->left, l, m, next);
            self(self, node->right, m, r, next);
        };
        visit(visit, root, 0, size_, id);
        return out.str();
    }
};
template <class S, class F> using PLazySegmentTree = PersistentLazySegmentTree<S, F>;

template <class V, class Op, class S, class Map, class Comp, class F>
auto initPersistentLazySegmentTree(const V &v, Op op, S e, Map map, Comp comp, F id) {
    return PersistentLazySegmentTree<S, F>(v, op, e, map, comp, id);
}

template <class V, class Op, class S, class Map, class Comp, class F>
auto initLazySegmentTree(const V &v, Op op, S e, Map map, Comp comp, F id) {
    return initPersistentLazySegmentTree(v, op, e, map, comp, id);
}

template <class... A> auto newPersistentLazySegWith(A &&...args) {
    return initPersistentLazySegmentTree(std::forward<A>(args)...);
}

template <class... A> auto newLazySegWith(A &&...args) {
    return initPersistentLazySegmentTree(std::forward<A>(args)...);
}

}
