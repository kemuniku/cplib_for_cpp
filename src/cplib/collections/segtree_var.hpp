#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <sstream>

namespace cplib {
// 構築O(N)、点更新・区間取得・境界探索O(log N)、全体取得O(1)。
template <class T, class Merge = std::function<T(T, T)>> class SegmentTree {
    T identity;
    [[no_unique_address]] Merge merge;
    Int lastnode = 1, length = 0;

public:
    struct Elem {
        SegmentTree *tree = nullptr;
        T v{};
        Int index = 0;

        // 対象の値を返す。O(1)。
        operator T() const {
            return v;
        }

        // 対象の値を更新する。O(log(n+1))回の演算。
        Elem &operator=(const T &value) {
            v = value;
            tree->propagate(index);
            return *this;
        }

        // 他の位置の値を読み取り代入する。O(log(n+1))回の演算。
        Elem &operator=(const Elem &value) {
            return *this = value.v;
        }

        // 未設定の要素参照を構築する。O(1)と値の初期化に要する時間。
        Elem() = default;

        // 木・値・木内位置を指定して要素参照を構築する。O(1)と値のコピーに要する時間。
        Elem(SegmentTree *t, T value, Int i) : tree(t), v(value), index(i) {
        }

        // 要素参照をコピーする。O(1)と値のコピーに要する時間。
        Elem(const Elem &) = default;
// 複合代入では値を更新して木に伝播する。O(log(n+1))回の演算と指定した値の演算。
#define CPLIB_SEG_ELEM_ASSIGN(OP)                                                                  \
    Elem &operator OP(const T &value) {                                                            \
        v OP value;                                                                                \
        tree->propagate(index);                                                                    \
        return *this;                                                                              \
    }

        CPLIB_SEG_ELEM_ASSIGN(+=)

        CPLIB_SEG_ELEM_ASSIGN(-=)
        CPLIB_SEG_ELEM_ASSIGN(*=)
        CPLIB_SEG_ELEM_ASSIGN(/=)
        CPLIB_SEG_ELEM_ASSIGN(^=)
        CPLIB_SEG_ELEM_ASSIGN(&=)
        CPLIB_SEG_ELEM_ASSIGN(|=)
        CPLIB_SEG_ELEM_ASSIGN(%=)
        CPLIB_SEG_ELEM_ASSIGN(>>=)
        CPLIB_SEG_ELEM_ASSIGN(<<=)
#undef CPLIB_SEG_ELEM_ASSIGN

        // 床除算して木に伝播する。O(log(n+1))回の木の演算と値の除算・剰余。
        // Nimの//=と**=に相当する名前付き演算。
        Elem &floor_div_assign(const T &value) {
            T q = v / value, r = v % value;
            if (r != 0 && ((r < 0) != (value < 0)))
                --q;
            v = q;
            tree->propagate(index);
            return *this;
        }

        // k乗して木に伝播する。O(log(k+1))回の乗算とO(log(n+1))回の木の演算。
        Elem &pow_assign(Int k) {
            T a = v, result = 1;
            assert(k >= 0);
            while (k) {
                if (k & 1)
                    result *= a;
                k >>= 1;
                if (k)
                    a *= a;
            }
            v = result;
            tree->propagate(index);
            return *this;
        }

        // 要素の値をストリームに出力する。計算量は値の出力に依存する。
        friend std::ostream &operator<<(std::ostream &s, const Elem &e) {
            return s << e.v;
        }
    };

    std::vector<Elem> arr;

    // 木内位置xから親へ積を再計算する。O(log(n+1))回の演算。
    void propagate(Int x) {
        while (x > 1) {
            x >>= 1;
            arr[x].v = merge(arr[2 * x].v, arr[2 * x + 1].v);
        }
    }

    // 長さnの単位元列から木を構築する。O(n+1)回の演算。
    SegmentTree(Int n, Merge op, T e) : identity(e), merge(op), length(n) {
        assert(n >= 0);
        while (lastnode < n)
            lastnode *= 2;
        arr.resize(2 * lastnode);
        for (Int i = 0; i < 2 * lastnode; ++i) {
            arr[i].tree = this;
            arr[i].index = i;
            arr[i].v = e;
        }
        for (Int i = lastnode - 1; i > 0; --i)
            arr[i].v = merge(arr[2 * i].v, arr[2 * i + 1].v);
    }

    // 配列から木を構築する。要素数をnとしてO(n+1)回の演算。
    SegmentTree(std::span<const T> v, Merge op, T e) : SegmentTree(v.size(), op, e) {
        for (Int i = 0; i < Int(v.size()); ++i)
            arr[lastnode + i].v = v[i];
        for (Int i = lastnode - 1; i > 0; --i)
            arr[i].v = merge(arr[2 * i].v, arr[2 * i + 1].v);
    }

    // 木をコピーする。要素数をnとしてO(n+1)。
    SegmentTree(const SegmentTree &) = default;
    // 木をムーブ構築する。O(1)。
    SegmentTree(SegmentTree &&) = default;
    // 木をムーブ代入する。代入先の要素数をnとしてO(n+1)。
    SegmentTree &operator=(SegmentTree &&) = default;

    // 木をコピーして代入する。両者の要素数の合計をnとしてO(n+1)。
    SegmentTree &operator=(const SegmentTree &other) {
        if (this != &other) {
            auto copy = other.arr;
            arr.swap(copy);
            identity = other.identity;
            merge = other.merge;
            lastnode = other.lastnode;
            length = other.length;
        }
        return *this;
    }

    // 要素数を返す。O(1)。
    Int len() const {
        return length;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return len();
    }

    // x番目の値を設定する。O(log(n+1))回の演算。
    void update(Int x, const T &value) {
        assert(0 <= x && x < length);
        x += lastnode;
        arr[x].v = value;
        propagate(x);
    }

    // O(log(n+1))回の演算。
    // 半開区間[l,r)の演算結果を返します。
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

    // 両端を含む区間の積を返す。O(log(n+1))回の演算。
    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return get(resolve_index(length, s.a), resolve_index(length, s.b) + 1);
    }

    // 両端を含む区間の積を返す。O(log(n+1))回の演算。
    template <class L, class R> T operator[](ClosedSlice<L, R> s) const {
        return get(s);
    }

    // i番目の値を返す。O(1)。
    T operator[](Int i) const {
        assert(0 <= i && i < length);
        return arr[i + lastnode];
    }

    // i番目を読み書きする参照を返す。O(1)。
    Elem &operator[](Int i) {
        assert(0 <= i && i < length);
        arr[lastnode + i].tree = this;
        return arr[lastnode + i];
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
    // 全区間の積を返す。O(1)。
    T get_all() const {
        return arr[1];
    }

    // 全要素を空白区切りで文字列化する。O(n+出力長)と各要素の文字列化の計算量。
    std::string str() const {
        std::ostringstream s;
        for (Int i = 0; i < length; ++i) {
            if (i)
                s << ' ';
            s << arr[lastnode + i];
        }
        return s.str();
    }

    // f(積[l,r))が真となる最大のrを返す。fは単位元で真かつ単調。O(log(n+1))回の演算と述語呼出し。
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

    // f(積[l,r))が真となる最小のlを返す。fは単位元で真かつ単調。O(log(n+1))回の演算と述語呼出し。
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

// 長さnの単位元列から木を構築する。O(n+1)回の演算。
template <class T, class Op> auto initSegmentTree(Int n, Op op, T e) {
    return SegmentTree<T>(n, op, e);
}

// 配列から木を構築する。要素数をnとしてO(n+1)回の演算。
template <class T, class Op> auto initSegmentTree(std::span<const T> v, Op op, T e) {
    return SegmentTree<T>(v, op, e);
}

// 配列から木を構築する。要素数をnとしてO(n+1)回の演算。
template <class T, class Op> auto initSegmentTree(const std::vector<T> &v, Op op, T e) {
    return SegmentTree<T>(v, op, e);
}

// 配列から木を構築する。要素数をnとしてO(n+1)回の演算。
template <class V, class Op, class T> auto newSegWith(const V &v, Op op, T e) {
    return cplib::initSegmentTree(v, op, e);
}

// 要素数を返す。O(1)。
template <class T, class M> Int len(const SegmentTree<T, M> &s) {
    return s.len();
}

// i番目の値を更新する。O(log(n+1))回の演算。
template <class T, class M> void update(SegmentTree<T, M> &s, Int i, const T &v) {
    s.update(i, v);
}

// 半開区間[l,r)の積を返す。O(log(n+1))回の演算。
template <class T, class M> auto get(const SegmentTree<T, M> &s, Int l, Int r) {
    return s.get(l, r);
}

// 両端を含む区間の積を返す。O(log(n+1))回の演算。
template <class T, class M, class L, class R>
auto get(const SegmentTree<T, M> &s, ClosedSlice<L, R> range) {
    return s.get(range);
}

// 全区間の積を返す。O(1)。
template <class T, class M> auto get_all(const SegmentTree<T, M> &s) {
    return s.get_all();
}

// 全要素を空白区切りで文字列化する。O(n+出力長)と各要素の文字列化の計算量。
template <class T, class M> auto to_string(const SegmentTree<T, M> &s) {
    return s.str();
}

// f(積[l,r))が真となる最大のrを返す。単調なfについてO(log(n+1))回の演算と述語呼出し。
template <class T, class M, class F> Int max_right(const SegmentTree<T, M> &s, Int l, F f) {
    return s.max_right(l, f);
}

// f(積[l,r))が真となる最小のlを返す。単調なfについてO(log(n+1))回の演算と述語呼出し。
template <class T, class M, class F> Int min_left(const SegmentTree<T, M> &s, Int r, F f) {
    return s.min_left(r, f);
}
}
