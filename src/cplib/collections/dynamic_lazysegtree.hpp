#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <memory>

namespace cplib {
// 区間をノードに持つ遅延AVL木。構築O(1)、各操作O(log(K+2))、Q回更新後O(Q+1)領域。
template <class S, class F> class DynamicLazySegmentTree {
    struct Node;
    using Ptr = std::shared_ptr<Node>;

    struct Node {
        Int a, b, lo, hi, height = 1;
        S value, product;
        F tag, lazy;
        bool pending = false;
        Ptr left, right;
    };

    Ptr root;
    Int length = 0, nodes = 0;
    std::function<S(S, S)> merge;
    S identity;
    std::function<S(F, S)> mapping;
    std::function<F(F, F)> composition;
    F id;
    std::function<S(Int, Int)> initial;

    Ptr makeNode(Int a, Int b, F tag) {
        S value = mapping(tag, initial(a, b));
        ++nodes;
        return std::make_shared<Node>(
            Node{a, b, a, b, 1, value, value, tag, id, false, nullptr, nullptr});
    }

    static Int height(const Ptr &p) {
        return p ? p->height : 0;
    }

    static Ptr clone(const Ptr &p) {
        if (!p)
            return nullptr;
        auto q = std::make_shared<Node>(*p);
        q->left = clone(p->left);
        q->right = clone(p->right);
        return q;
    }

    void pull(const Ptr &p) {
        p->product = p->value;
        p->lo = p->a;
        p->hi = p->b;
        if (p->left) {
            p->product = merge(p->left->product, p->product);
            p->lo = p->left->lo;
        }
        if (p->right) {
            p->product = merge(p->product, p->right->product);
            p->hi = p->right->hi;
        }
        p->height = std::max(height(p->left), height(p->right)) + 1;
    }

    void allApply(const Ptr &p, F f) {
        if (!p)
            return;
        p->value = mapping(f, p->value);
        p->product = mapping(f, p->product);
        p->tag = composition(f, p->tag);
        if (p->pending)
            p->lazy = composition(f, p->lazy);
        else {
            p->lazy = f;
            p->pending = true;
        }
    }

    void push(const Ptr &p) {
        if (p->pending) {
            allApply(p->left, p->lazy);
            allApply(p->right, p->lazy);
            p->lazy = id;
            p->pending = false;
        }
    }

    Ptr rotateLeft(Ptr p) {
        push(p);
        auto result = p->right;
        push(result);
        p->right = result->left;
        result->left = p;
        pull(p);
        pull(result);
        return result;
    }

    Ptr rotateRight(Ptr p) {
        push(p);
        auto result = p->left;
        push(result);
        p->left = result->right;
        result->right = p;
        pull(p);
        pull(result);
        return result;
    }

    Ptr balance(Ptr p) {
        pull(p);
        if (height(p->left) > height(p->right) + 1) {
            if (height(p->left->left) < height(p->left->right))
                p->left = rotateLeft(p->left);
            return rotateRight(p);
        }
        if (height(p->right) > height(p->left) + 1) {
            if (height(p->right->right) < height(p->right->left))
                p->right = rotateRight(p->right);
            return rotateLeft(p);
        }
        return p;
    }

    Ptr insertFirst(Ptr p, Ptr added) {
        if (!p)
            return added;
        push(p);
        p->left = insertFirst(p->left, added);
        return balance(p);
    }

    Ptr splitAt(Ptr p, Int index) {
        if (!p || index == p->a || index == p->b)
            return p;
        push(p);
        if (index < p->a)
            p->left = splitAt(p->left, index);
        else if (p->b < index)
            p->right = splitAt(p->right, index);
        else {
            auto added = makeNode(index, p->b, p->tag);
            p->b = index;
            p->value = mapping(p->tag, initial(p->a, index));
            p->right = insertFirst(p->right, added);
        }
        return balance(p);
    }

    void applyNode(const Ptr &p, Int l, Int r, F f) {
        if (!p || r <= p->lo || p->hi <= l)
            return;
        if (l <= p->lo && p->hi <= r) {
            allApply(p, f);
            return;
        }
        push(p);
        applyNode(p->left, l, r, f);
        if (l <= p->a && p->b <= r) {
            p->value = mapping(f, p->value);
            p->tag = composition(f, p->tag);
        }
        applyNode(p->right, l, r, f);
        pull(p);
    }

    S getNode(const Ptr &p, Int l, Int r) {
        if (!p || r <= p->lo || p->hi <= l)
            return identity;
        if (l <= p->lo && p->hi <= r)
            return p->product;
        push(p);
        S result = getNode(p->left, l, r);
        Int a = std::max(l, p->a), b = std::min(r, p->b);
        if (a < b) {
            S value = a == p->a && b == p->b ? p->value : mapping(p->tag, initial(a, b));
            result = merge(result, value);
        }
        return merge(result, getNode(p->right, l, r));
    }

    void setNode(const Ptr &p, Int index, S value) {
        push(p);
        if (index < p->a)
            setNode(p->left, index, value);
        else if (index >= p->b)
            setNode(p->right, index, value);
        else {
            p->value = value;
            p->tag = id;
        }
        pull(p);
    }

public:
    // [0,n)をO(1)で生成します。comp(f,g)はgの後にfを適用します。
    // init(l,r)は更新を一度も行っていない配列の半開区間[l,r)の集約値です。
    // 例えばS=(sum,size)の零初期化では(0,r-l)を返し、空区間の単位元(0,0)とは区別します。
    // 区間分割や部分取得にも使い、保存済みの作用は木が適用するため、initに更新を反映させないでください。
    // init(l,r)=op(init(l,m),init(m,r))、init(l,l)=eを満たす必要があります。
    // 記載の計算量はinitを含む各コールバックがO(1)の場合です。
    template <class Op, class Map, class Comp, class Initial>
    DynamicLazySegmentTree(Int n, Op op, S e, Map map, Comp comp, F id_, Initial init)
        : length(n), merge(op), identity(e), mapping(map), composition(comp), id(id_),
          initial(init) {
        assert(n >= 0);
        if (n > 0)
            root = makeNode(0, n, id);
    }

    // 木を複製する。保持区間数をKとしてO(K)。
    DynamicLazySegmentTree(const DynamicLazySegmentTree &p)
        : root(clone(p.root)), length(p.length), nodes(p.nodes), merge(p.merge),
          identity(p.identity), mapping(p.mapping), composition(p.composition), id(p.id),
          initial(p.initial) {
    }

    // 木をムーブ構築する。O(1)。
    DynamicLazySegmentTree(DynamicLazySegmentTree &&) = default;
    // 木をムーブ代入する。代入先の保持区間数をKとしてO(K)。
    DynamicLazySegmentTree &operator=(DynamicLazySegmentTree &&) = default;

    // 木を複製して代入する。両者の保持区間数の合計をKとしてO(K)。
    DynamicLazySegmentTree &operator=(const DynamicLazySegmentTree &p) {
        if (this != &p) {
            DynamicLazySegmentTree copy(p);
            *this = std::move(copy);
        }
        return *this;
    }

    // 更新は両端の境界を作る。追加ノードは高々2個。
    // [l,r)へ最悪O(log(K+2))で作用させます。追加ノードは高々2個、Q回更新後の空間はO(Q+1)。
    void apply(Int l, Int r, F f) {
        assert(0 <= l && l <= r && r <= length);
        if (l == r)
            return;
        root = splitAt(root, l);
        root = splitAt(root, r);
        applyNode(root, l, r, f);
    }

    // 半開区間[l,r)の積を最悪O(log(K+2))で返します。取得ではノードを追加しません。
    S get(Int l, Int r) {
        assert(0 <= l && l <= r && r <= length);
        return l == r ? identity : getNode(root, l, r);
    }

    // p番目の値を更新する。保持区間数をKとしてO(log(K+2))。
    void update(Int p, S value) {
        assert(0 <= p && p < length);
        root = splitAt(root, p);
        root = splitAt(root, p + 1);
        setNode(root, p, value);
    }

    // 両端を含む区間の積を返す。保持区間数をKとしてO(log(K+2))。
    template <class L, class R> S get(ClosedSlice<L, R> s) {
        return get(resolve_index(length, s.a), resolve_index(length, s.b) + 1);
    }

    // 両端を含む区間にfを作用させる。保持区間数をKとしてO(log(K+2))。
    template <class L, class R> void apply(ClosedSlice<L, R> s, F f) {
        apply(resolve_index(length, s.a), resolve_index(length, s.b) + 1, f);
    }

    // 両端を含む区間の積を返す。保持区間数をKとしてO(log(K+2))。
    template <class L, class R> S operator[](ClosedSlice<L, R> s) {
        return get(s);
    }

    struct Reference {
        DynamicLazySegmentTree *tree;
        Int index;

        // 対象の値を取得する。保持区間数をKとしてO(log(K+2))。
        operator S() const {
            return tree->get(index, index + 1);
        }

        // 対象の値を設定する。保持区間数をKとしてO(log(K+2))。
        Reference &operator=(S value) {
            tree->update(index, value);
            return *this;
        }

        // 他の位置の値を読み取り代入する。保持区間数をKとしてO(log(K+2))。
        Reference &operator=(const Reference &r) {
            return *this = S(r);
        }
    };

    // p番目を読み書きする参照を返す。O(1)。
    Reference operator[](Int p) {
        assert(0 <= p && p < length);
        return {this, p};
    }

    // 末尾からp番目を読み書きする参照を返す。O(1)。
    Reference operator[](BackwardsIndex p) {
        return (*this)[length - p.value];
    }

    // 要素数を返す。O(1)。
    Int len() const {
        return length;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return len();
    }

    // 保持する区間数K（確保したノード数）をO(1)で返します。
    Int node_count() const {
        return nodes;
    }

    // 全区間の積を返す。O(1)。
    S get_all() const {
        return root ? root->product : identity;
    }
};

// 長さと演算・作用・初期区間の積を指定して構築する。各コールバックがO(1)ならO(1)。
template <class S, class F, class Op, class Map, class Comp, class Initial>
auto initDynamicLazySegmentTree(Int n, Op op, S e, Map map, Comp comp, F id, Initial initial) {
    return DynamicLazySegmentTree<S, F>(n, op, e, map, comp, id, initial);
}

// 長さと演算・作用・初期区間の積を指定して構築する。各コールバックがO(1)ならO(1)。
template <class... Args> auto newDynamicLazySegWith(Args &&...args) {
    return initDynamicLazySegmentTree(std::forward<Args>(args)...);
}

// 要素数を返す。O(1)。
template <class S, class F> Int len(const DynamicLazySegmentTree<S, F> &t) {
    return t.len();
}

// 保持する区間数を返す。O(1)。
template <class S, class F> Int node_count(const DynamicLazySegmentTree<S, F> &t) {
    return t.node_count();
}

// 全区間の積を返す。O(1)。
template <class S, class F> S get_all(const DynamicLazySegmentTree<S, F> &t) {
    return t.get_all();
}

// 半開区間[l,r)の積を返す。保持区間数をKとしてO(log(K+2))。
template <class S, class F> S get(DynamicLazySegmentTree<S, F> &t, Int l, Int r) {
    return t.get(l, r);
}

// 両端を含む区間の積を返す。保持区間数をKとしてO(log(K+2))。
template <class S, class F, class L, class R>
S get(DynamicLazySegmentTree<S, F> &t, ClosedSlice<L, R> s) {
    return t.get(s);
}

// 半開区間[l,r)にfを作用させる。保持区間数をKとしてO(log(K+2))。
template <class S, class F> void apply(DynamicLazySegmentTree<S, F> &t, Int l, Int r, F f) {
    t.apply(l, r, f);
}

// 両端を含む区間にfを作用させる。保持区間数をKとしてO(log(K+2))。
template <class S, class F, class L, class R>
void apply(DynamicLazySegmentTree<S, F> &t, ClosedSlice<L, R> s, F f) {
    t.apply(s, f);
}

// p番目の値を設定する。保持区間数をKとしてO(log(K+2))。
template <class S, class F> void update(DynamicLazySegmentTree<S, F> &t, Int p, S value) {
    t.update(p, value);
}
}
