#pragma once
#include <cplib/collections/lazysegtree.hpp>
#include <cplib/collections/compressed_coordinates_internal.hpp>

namespace cplib {
// 登録点モードと隣接座標間の区間モード。構築O(N log N)、操作O(log N)。
template <class K, class S, class F> class CompressedLazySegmentTree {
    using Tree = LazySegmentTree<S, F>;
    std::vector<K> coords;
    Tree tree;
    S identity;
    F actionIdentity;
    bool intervals = false;
    std::vector<Int> indexSlots;
    std::function<S(Tree &, Int, Int)> getImpl;
    std::function<void(Tree &, Int, Int, F)> applyImpl;
    std::function<void(Tree &, Int, S)> updateImpl;
    std::function<S(Tree &, Int)> pointImpl;

    Int boundary(const K &x) const {
        if (intervals) {
            Int i = detail::findCompressedCoordinate<K>(coords, indexSlots, x);
            assert(i >= 0);
            return i;
        }
        return std::lower_bound(coords.begin(), coords.end(), x) - coords.begin();
    }

    std::pair<Int, Int> sliceBounds(ClosedSlice<K, K> s) const {
        if (s.b < s.a)
            return {0, 0};
        if (intervals) {
            if constexpr (std::integral<K>) {
                assert(s.b < std::numeric_limits<K>::max());
                return {boundary(s.a), boundary(K(s.b + 1))};
            } else {
                assert(false && "interval slices require ordinal coordinates");
                return {0, 0};
            }
        }
        return {std::lower_bound(coords.begin(), coords.end(), s.a) - coords.begin(),
                std::upper_bound(coords.begin(), coords.end(), s.b) - coords.begin()};
    }

    S range(Int a, Int b) {
        return getImpl ? getImpl(tree, a, b) : tree.get(a, b);
    }

    void rangeApply(Int a, Int b, F f) {
        if (applyImpl)
            applyImpl(tree, a, b, f);
        else
            tree.apply(a, b, f);
    }

public:
    // 登録座標は構築後に追加できません。Kには一貫した<と==が必要です。
    // 点モードでは各点をinitial(x)（省略時はe）で初期化します。
    // 区間モードでは各葉[x[i],x[i+1])をinitial(x[i],x[i+1])で初期化します。
    // 区間モードの操作端点は事前登録してください。区間和ではSに圧縮前の区間長を含めます。
    // comp(f,g)はgの後にfを適用します。initialを含むコールバックはO(1)を想定します。
    template <class Op, class Map, class Comp, class Initial = std::nullptr_t>
    CompressedLazySegmentTree(std::span<const K> input, Op op, S e, Map map, Comp comp, F id,
                              Initial initial = nullptr)
        : coords(detail::uniqueCompressedCoordinates(input)),
          tree(Int(coords.size()), op, e, map, comp, id), identity(e), actionIdentity(id),
          indexSlots(detail::initCompressedCoordinateIndex<K>(coords)) {
        if constexpr (!std::is_same_v<Initial, std::nullptr_t>) {
            std::vector<S> values;
            if constexpr (std::is_invocable_v<Initial, K, K>) {
                intervals = true;
                for (Int i = 0; i + 1 < Int(coords.size()); ++i)
                    values.push_back(initial(coords[i], coords[i + 1]));
            } else {
                for (const K &x : coords)
                    values.push_back(initial(x));
            }
            tree = initLazySegmentTree(values, op, e, map, comp, id);
        }
    }

    // 半開区間[l,r)の積をO(log N)で返します。登録点を保持する場合は両端が未登録でも構いません。
    S get(const K &l, const K &r) {
        assert(!(r < l));
        return range(boundary(l), boundary(r));
    }

    // 閉区間の積をO(log N)で返します。逆順の区間では単位元を返します。
    // 区間モードではKは順序型とし、左端と右端の次の座標を事前登録してください。
    S get(ClosedSlice<K, K> s) {
        auto [a, b] = sliceBounds(s);
        return range(a, b);
    }

    // 半開区間[l,r)へO(log N)で作用させます。登録点を保持する場合は未登録点には作用しません。
    void apply(const K &l, const K &r, F f) {
        assert(!(r < l));
        rangeApply(boundary(l), boundary(r), f);
    }

    // 閉区間へO(log N)で作用させます。逆順の区間では何もしません。
    // 区間モードではKは順序型とし、左端と右端の次の座標を事前登録してください。
    void apply(ClosedSlice<K, K> s, F f) {
        auto [a, b] = sliceBounds(s);
        rangeApply(a, b, f);
    }

    // 登録点xの値をO(log N)で上書きします。区間を葉とする場合は[x,次の登録座標)全体を上書きします。
    void update(const K &x, S value) {
        Int i = detail::findCompressedCoordinate<K>(coords, indexSlots, x);
        assert(i >= 0 && (!intervals || i + 1 < Int(coords.size())));
        if (updateImpl)
            updateImpl(tree, i, value);
        else
            tree.update(i, value);
    }

    // 登録点xの値をO(log N)で返します。区間を葉とする場合は[x,次の登録座標)の積を返します。
    // 登録点を保持する場合、未登録点の値は単位元です。区間を葉とする場合は未登録座標を許しません。
    S point(const K &x) {
        Int i = detail::findCompressedCoordinate<K>(coords, indexSlots, x);
        if (intervals)
            assert(i >= 0 && i + 1 < Int(coords.size()));
        else if (i < 0)
            return identity;
        return pointImpl ? pointImpl(tree, i) : tree.get(i);
    }

    struct Reference {
        CompressedLazySegmentTree *tree;
        K key;

        operator S() const {
            return tree->point(key);
        }

        Reference &operator=(const S &v) {
            tree->update(key, v);
            return *this;
        }

        Reference &operator=(const Reference &v) {
            return *this = S(v);
        }
    };

    Reference operator[](const K &x) {
        return {this, x};
    }

    S operator[](ClosedSlice<K, K> s) {
        return get(s);
    }

    S get_all() const {
        return tree.get_all();
    }

    Int len() const {
        return intervals ? std::max<Int>(0, Int(coords.size()) - 1) : Int(coords.size());
    }

    Int size() const {
        return len();
    }

    std::string str() {
        return tree.str();
    }

    std::pair<S, F> compressedLazyDefaults() const {
        return {identity, actionIdentity};
    }

    template <class Get, class Apply, class Update, class Point>
    void setCompressedLazyOperations(Get get, Apply apply, Update update, Point point) {
        getImpl = get;
        applyImpl = apply;
        updateImpl = update;
        pointImpl = point;
    }
};

template <class K, class S, class F, class Op, class Map, class Comp,
          class Initial = std::nullptr_t>
auto initCompressedLazySegmentTree(std::span<const K> coords, Op op, S e, Map map, Comp comp, F id,
                                   Initial initial = nullptr) {
    return CompressedLazySegmentTree<K, S, F>(coords, op, e, map, comp, id, initial);
}

template <class K, class S, class F, class Op, class Map, class Comp,
          class Initial = std::nullptr_t>
auto initCompressedLazySegmentTree(const std::vector<K> &coords, Op op, S e, Map map, Comp comp,
                                   F id, Initial initial = nullptr) {
    return CompressedLazySegmentTree<K, S, F>(coords, op, e, map, comp, id, initial);
}

// 直接演算版。元のマクロと同じ境界の共有祖先省略・関数特殊化を保持する。
template <class Coords, class Op, class S, class Map, class Comp, class F,
          class Initial = std::nullptr_t>
auto newCompressedLazySegWith(const Coords &coords, Op op, S e, Map map, Comp comp, F id,
                              Initial initial = nullptr) {
    using Tree = LazySegmentTree<S, F>;
    auto result = initCompressedLazySegmentTree(coords, op, e, map, comp, id, initial);
    auto applyNode = [map, comp](Tree &t, Int node, F f) {
        t.arr[node] = map(f, t.arr[node]);
        if (node < Int(t.lazy.size()))
            t.lazy[node] = comp(f, t.lazy[node]);
    };
    auto pushNode = [applyNode, id](Tree &t, Int node) {
        applyNode(t, node * 2, t.lazy[node]);
        applyNode(t, node * 2 + 1, t.lazy[node]);
        t.lazy[node] = id;
    };
    auto pushBoundary = [pushNode](Tree &t, Int l, Int r) {
        int lz = std::countr_zero(UInt(l)), rz = std::countr_zero(UInt(r));
        for (int i = std::countr_zero(t.lazy.size()); i >= 1; --i) {
            Int lp = l >> i, rp = (r - 1) >> i;
            if (i > lz)
                pushNode(t, lp);
            if (i > rz && (i <= lz || lp != rp))
                pushNode(t, rp);
        }
    };
    auto get = [pushBoundary, op, e](Tree &t, Int a, Int b) {
        if (a == b)
            return e;
        Int l = a + t.lazy.size(), r = b + t.lazy.size();
        pushBoundary(t, l, r);
        S left = e, right = e;
        while (l < r) {
            if (l & 1)
                left = op(left, t.arr[l++]);
            if (r & 1)
                right = op(t.arr[--r], right);
            l >>= 1;
            r >>= 1;
        }
        return op(left, right);
    };
    auto apply = [pushBoundary, applyNode, op](Tree &t, Int a, Int b, F f) {
        if (a == b)
            return;
        Int first = a + t.lazy.size(), last = b + t.lazy.size();
        pushBoundary(t, first, last);
        Int l = first, r = last;
        while (l < r) {
            if (l & 1)
                applyNode(t, l++, f);
            if (r & 1)
                applyNode(t, --r, f);
            l >>= 1;
            r >>= 1;
        }
        int lz = std::countr_zero(UInt(first)), rz = std::countr_zero(UInt(last));
        for (int i = 1; i <= std::countr_zero(t.lazy.size()); ++i) {
            Int lp = first >> i, rp = (last - 1) >> i;
            if (i > lz)
                t.arr[lp] = op(t.arr[lp * 2], t.arr[lp * 2 + 1]);
            if (i > rz && (i <= lz || lp != rp))
                t.arr[rp] = op(t.arr[rp * 2], t.arr[rp * 2 + 1]);
        }
    };
    auto update = [pushNode, op](Tree &t, Int i, S value) {
        Int node = i + t.lazy.size();
        for (int level = std::countr_zero(t.lazy.size()); level >= 1; --level)
            pushNode(t, node >> level);
        t.arr[node] = value;
        while (node > 1) {
            node >>= 1;
            t.arr[node] = op(t.arr[2 * node], t.arr[2 * node + 1]);
        }
    };
    auto point = [pushNode](Tree &t, Int i) {
        Int node = i + t.lazy.size();
        for (int level = std::countr_zero(t.lazy.size()); level >= 1; --level)
            pushNode(t, node >> level);
        return t.arr[node];
    };
    result.setCompressedLazyOperations(get, apply, update, point);
    return result;
}

template <class K, class S, class F> Int len(const CompressedLazySegmentTree<K, S, F> &t) {
    return t.len();
}

template <class K, class S, class F>
S get(CompressedLazySegmentTree<K, S, F> &t, const K &l, const K &r) {
    return t.get(l, r);
}

template <class K, class S, class F>
S get(CompressedLazySegmentTree<K, S, F> &t, ClosedSlice<K, K> r) {
    return t.get(r);
}

template <class K, class S, class F> S get_all(const CompressedLazySegmentTree<K, S, F> &t) {
    return t.get_all();
}

template <class K, class S, class F>
void apply(CompressedLazySegmentTree<K, S, F> &t, const K &l, const K &r, F f) {
    t.apply(l, r, f);
}

template <class K, class S, class F>
void apply(CompressedLazySegmentTree<K, S, F> &t, ClosedSlice<K, K> r, F f) {
    t.apply(r, f);
}

template <class K, class S, class F>
void update(CompressedLazySegmentTree<K, S, F> &t, const K &x, S value) {
    t.update(x, value);
}

template <class K, class S, class F> auto to_string(CompressedLazySegmentTree<K, S, F> &t) {
    return t.str();
}
}
