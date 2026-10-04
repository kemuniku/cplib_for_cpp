#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <sstream>

namespace cplib::detail::segtree_base {
// 構築O(N)、点更新・区間取得・境界探索O(log N)、全体取得O(1)。
template <class T, class Merge = std::function<T(T, T)>> class SegmentTree {
    T identity;
    [[no_unique_address]] Merge merge;
    Int lastnode = 1, length = 0;

public:
    std::vector<T> arr;

    // 長さnの全要素を単位元eで構築する。opに区間のマージ関数を指定する。
    SegmentTree(Int n, Merge op, T e) : identity(e), merge(op), length(n) {
        assert(n >= 0);
        while (lastnode < n)
            lastnode *= 2;
        arr.assign(2 * lastnode, e);
        for (Int i = lastnode - 1; i > 0; --i)
            arr[i] = merge(arr[2 * i], arr[2 * i + 1]);
    }

    // vを元に構築する。opに区間のマージ関数、eに単位元を指定する。
    SegmentTree(std::span<const T> v, Merge op, T e) : SegmentTree(v.size(), op, e) {
        std::copy(v.begin(), v.end(), arr.begin() + lastnode);
        for (Int i = lastnode - 1; i > 0; --i)
            arr[i] = merge(arr[2 * i], arr[2 * i + 1]);
    }

    Int len() const {
        return length;
    }

    Int size() const {
        return len();
    }

    // xの要素をvalueに変更する。
    void update(Int x, const T &value) {
        assert(0 <= x && x < length);
        x += lastnode;
        arr[x] = value;
        while (x > 1) {
            x >>= 1;
            arr[x] = merge(arr[2 * x], arr[2 * x + 1]);
        }
    }

    // 半開区間[l,r)の演算結果を返す。
    T get(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= length);
        l += lastnode;
        r += lastnode;
        T a = identity, b = identity;
        while (l < r) {
            if (l & 1)
                a = merge(a, arr[l++]);
            if (r & 1)
                b = merge(arr[--r], b);
            l >>= 1;
            r >>= 1;
        }
        return merge(a, b);
    }

    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return get(resolve_index(length, s.a), resolve_index(length, s.b) + 1);
    }

    template <class L, class R> T operator[](ClosedSlice<L, R> s) const {
        return get(s);
    }

    T operator[](Int i) const {
        assert(0 <= i && i < length);
        return arr[i + lastnode];
    }

    struct Reference {
        SegmentTree *tree;
        Int index;

        operator T() const {
            return std::as_const(*tree)[index];
        }

        Reference &operator=(const T &v) {
            tree->update(index, v);
            return *this;
        }

        Reference &operator=(const Reference &v) {
            return *this = T(v);
        }
    };

    Reference operator[](Int i) {
        assert(0 <= i && i < length);
        return {this, i};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS

    // 全区間[0,len())の演算結果をO(1)で返す。
    T get_all() const {
        return arr[1];
    }

    std::string str() const {
        std::ostringstream s;
        for (Int i = 0; i < length; ++i) {
            if (i)
                s << ' ';
            s << arr[lastnode + i];
        }
        return s.str();
    }

    template <class F> Int max_right(Int l, F f) const {
        assert(0 <= l && l <= length && f(identity));
        if (l == length)
            return length;
        l += lastnode;
        T sm = identity;
        do {
            while (l % 2 == 0)
                l >>= 1;
            if (!f(merge(sm, arr[l]))) {
                while (l < lastnode) {
                    l *= 2;
                    if (f(merge(sm, arr[l]))) {
                        sm = merge(sm, arr[l]);
                        ++l;
                    }
                }
                return l - lastnode;
            }
            sm = merge(sm, arr[l]);
            ++l;
        } while ((l & -l) != l);
        return length;
    }

    template <class F> Int min_left(Int r, F f) const {
        assert(0 <= r && r <= length && f(identity));
        if (r == 0)
            return 0;
        r += lastnode;
        T sm = identity;
        do {
            --r;
            while (r > 1 && r % 2 != 0)
                r >>= 1;
            if (!f(merge(arr[r], sm))) {
                while (r < lastnode) {
                    r = 2 * r + 1;
                    if (f(merge(arr[r], sm))) {
                        sm = merge(arr[r], sm);
                        --r;
                    }
                }
                return r + 1 - lastnode;
            }
            sm = merge(arr[r], sm);
        } while ((r & -r) != r);
        return 0;
    }
};

template <class T, class Op> auto initSegmentTree(Int n, Op op, T e) {
    return SegmentTree<T>(n, op, e);
}

template <class T, class Op> auto initSegmentTree(std::span<const T> v, Op op, T e) {
    return SegmentTree<T>(v, op, e);
}

template <class T, class Op> auto initSegmentTree(const std::vector<T> &v, Op op, T e) {
    return SegmentTree<T>(v, op, e);
}

template <class V, class Op, class T> auto newSegWith(const V &v, Op op, T e) {
    return initSegmentTree(v, op, e);
}

template <class T, class M> Int len(const SegmentTree<T, M> &s) {
    return s.len();
}

template <class T, class M> void update(SegmentTree<T, M> &s, Int i, const T &v) {
    s.update(i, v);
}

template <class T, class M> auto get(const SegmentTree<T, M> &s, Int l, Int r) {
    return s.get(l, r);
}

template <class T, class M, class L, class R>
auto get(const SegmentTree<T, M> &s, ClosedSlice<L, R> range) {
    return s.get(range);
}

template <class T, class M> auto get_all(const SegmentTree<T, M> &s) {
    return s.get_all();
}

template <class T, class M> auto to_string(const SegmentTree<T, M> &s) {
    return s.str();
}

template <class T, class M, class F> Int max_right(const SegmentTree<T, M> &s, Int l, F f) {
    return s.max_right(l, f);
}

template <class T, class M, class F> Int min_left(const SegmentTree<T, M> &s, Int r, F f) {
    return s.min_left(r, f);
}
}
