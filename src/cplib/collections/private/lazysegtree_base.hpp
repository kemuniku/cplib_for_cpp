#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <bit>
#include <sstream>

namespace cplib::detail::lazysegtree_base {
// 通常版は単位作用、静的版はhasLazyフラグを使う別の特殊化。構築O(N)、操作O(log N)。
template <class S, class F, class Merge = std::function<S(S, S)>,
          class Mapping = std::function<S(F, S)>, class Composition = std::function<F(F, F)>,
          bool HasFlag = false>
class LazySegmentTree {
    S identity;
    F id;
    [[no_unique_address]] Merge merge;
    [[no_unique_address]] Mapping mapping;
    [[no_unique_address]] Composition composition;
    Int lastnode = 1, log = 0, length = 0;
    std::vector<bool> hasLazy;

    void pull(Int p) {
        arr[p] = merge(arr[2 * p], arr[2 * p + 1]);
    }

    void all_apply(Int p, const F &f) {
        arr[p] = mapping(f, arr[p]);
        if (p < lastnode) {
            if constexpr (HasFlag) {
                if (hasLazy[p])
                    lazy[p] = composition(f, lazy[p]);
                else {
                    lazy[p] = f;
                    hasLazy[p] = true;
                }
            } else
                lazy[p] = composition(f, lazy[p]);
        }
    }

    [[gnu::noinline]] void pushFlag(Int p) {
        if (hasLazy[p]) {
            F f = lazy[p];
            all_apply(2 * p, f);
            all_apply(2 * p + 1, f);
            hasLazy[p] = false;
        }
    }

    void push(Int p) {
        if constexpr (HasFlag)
            pushFlag(p);
        else {
            all_apply(2 * p, lazy[p]);
            all_apply(2 * p + 1, lazy[p]);
            lazy[p] = id;
        }
    }

    void all_push(Int p) {
        for (Int i = log; i >= 1; --i)
            push(p >> i);
    }

    // 静的版は左右が共有する祖先を一度だけ処理する。
    [[gnu::noinline]] void pushFlagBoundaries(Int l, Int r) {
        int lz = std::countr_zero(UInt(l)), rz = std::countr_zero(UInt(r));
        for (Int i = log; i >= 1; --i) {
            Int left = 0;
            if (i > lz) {
                left = l >> i;
                push(left);
            }
            if (i > rz) {
                Int right = (r - 1) >> i;
                if (right != left)
                    push(right);
            }
        }
    }

    void pushBoundaries(Int l, Int r) {
        if constexpr (HasFlag)
            pushFlagBoundaries(l, r);
        else {
            for (Int i = log; i > std::countr_zero(UInt(l)); --i)
                push(l >> i);
            for (Int i = log; i > std::countr_zero(UInt(r)); --i)
                push((r - 1) >> i);
        }
    }

    [[gnu::noinline]] void pullFlagBoundaries(Int l, Int r) {
        Int lm = std::countr_zero(UInt(l)) + 1, rm = std::countr_zero(UInt(r)) + 1;
        for (Int i = std::min(lm, rm); i <= log; ++i) {
            Int left = 0;
            if (i >= lm) {
                left = l >> i;
                pull(left);
            }
            if (i >= rm) {
                Int right = (r - 1) >> i;
                if (right != left)
                    pull(right);
            }
        }
    }

    void pullBoundaries(Int l, Int r) {
        if constexpr (HasFlag)
            pullFlagBoundaries(l, r);
        else {
            for (Int i = std::countr_zero(UInt(l)) + 1; i <= log; ++i)
                pull(l >> i);
            for (Int i = std::countr_zero(UInt(r)) + 1; i <= log; ++i)
                pull((r - 1) >> i);
        }
    }

public:
    using value_type = S;
    using action_type = F;

    // モノイドの単位元を返す。O(1)と値のコピーに要する時間。
    S calc_e() const {
        return identity;
    }

    // 作用の単位元を返す。O(1)と値のコピーに要する時間。
    F calc_id() const {
        return id;
    }

    // 指定したモノイド演算を呼び出す。計算量はその演算に依存する。
    S calc_op(const S &a, const S &b) const {
        return merge(a, b);
    }

    // 指定した作用の合成を呼び出す。計算量は合成に依存する。
    F calc_composition(const F &a, const F &b) const {
        return composition(a, b);
    }

    // 全葉への遅延伝播はO(N)。functional graphの全頂点集約に用いる。
    std::vector<S> materialize() {
        for (Int p = 1; p < lastnode; ++p)
            push(p);
        return std::vector<S>(arr.begin() + lastnode, arr.begin() + lastnode + length);
    }

    // [l,r)の葉をoutputに追加する。O(r-l+log(n+1))回の作用・合成。
    // 区間の葉を一括走査し、点取得の繰り返しによるlog N倍を避ける。
    void append_range(Int l, Int r, std::vector<S> &output) {
        assert(0 <= l && l <= r && r <= length);
        auto visit = [&](auto &&self, Int p, Int a, Int b) -> void {
            if (r <= a || b <= l)
                return;
            if (b - a == 1) {
                output.push_back(arr[p]);
                return;
            }
            push(p);
            Int mid = (a + b) / 2;
            self(self, 2 * p, a, mid);
            self(self, 2 * p + 1, mid, b);
        };
        if (l < r)
            visit(visit, 1, 0, lastnode);
    }

    std::vector<S> arr;
    std::vector<F> lazy;

    // 長さnの単位元列から構築する。O(n+1)時間・領域。
    LazySegmentTree(Int n, Merge op, S e, Mapping map, Composition comp, F id_)
        : identity(e), id(id_), merge(op), mapping(map), composition(comp), length(n) {
        assert(n >= 0);
        while (lastnode < n) {
            lastnode *= 2;
            ++log;
        }
        arr.assign(2 * lastnode, e);
        if constexpr (HasFlag) {
            lazy.resize(lastnode);
            hasLazy.assign(lastnode, false);
        } else
            lazy.assign(lastnode, id);
    }

    // 配列から構築する。要素数をnとしてO(n+1)回の演算。
    LazySegmentTree(std::span<const S> v, Merge op, S e, Mapping map, Composition comp, F id_)
        : LazySegmentTree(v.size(), op, e, map, comp, id_) {
        std::copy(v.begin(), v.end(), arr.begin() + lastnode);
        for (Int i = lastnode - 1; i > 0; --i)
            pull(i);
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return length;
    }

    // p番目の値を更新する。O(log(n+1))回の演算・作用・合成。
    void update(Int p, const S &value) {
        assert(0 <= p && p < length);
        p += lastnode;
        all_push(p);
        arr[p] = value;
        for (Int i = 1; i <= log; ++i)
            pull(p >> i);
    }

    // p番目の値を返す。O(log(n+1))回の作用・合成。
    S get(Int p) {
        assert(0 <= p && p < length);
        all_push(p + lastnode);
        return arr[p + lastnode];
    }

    // O(log(n+1))回の演算・作用・合成。
    // 半開区間[l,r)の積を返します。
    S get(Int l, Int r) {
        assert(0 <= l && l <= r && r <= length);
        if (l == r)
            return identity;
        l += lastnode;
        r += lastnode;
        pushBoundaries(l, r);
        S a = identity, b = identity;
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

    // 全要素の積を返す。O(1)。
    S get_all() const {
        return arr[1];
    }

    // 両端を含む区間の積を返す。O(log(n+1))回の演算・作用・合成。
    template <class L, class R> S get(ClosedSlice<L, R> s) {
        return get(resolve_index(length, s.a), resolve_index(length, s.b) + 1);
    }

    // 両端を含む区間の積を返す。O(log(n+1))回の演算・作用・合成。
    template <class L, class R> S operator[](ClosedSlice<L, R> s) {
        return get(s);
    }

    struct Reference {
        LazySegmentTree *tree;
        Int index;

        // 対象の値を取得する。O(log(n+1))回の作用・合成。
        operator S() const {
            return tree->get(index);
        }

        // 対象の値を更新する。O(log(n+1))回の演算・作用・合成。
        Reference &operator=(const S &v) {
            tree->update(index, v);
            return *this;
        }

        // 他の位置の値を読み取り代入する。O(log(n+1))回の演算・作用・合成。
        Reference &operator=(const Reference &v) {
            return *this = S(v);
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

    // 全要素を空白区切りで文字列化する。O(n log(n+1))回の作用・合成と出力処理。
    std::string str() {
        std::ostringstream s;
        for (Int i = 0; i < length; ++i) {
            if (i)
                s << ' ';
            s << get(i);
        }
        return s.str();
    }

    // O(log(n+1))回の演算・作用・合成。
    // 半開区間[l,r)にfを作用させます。
    void apply(Int l, Int r, const F &f) {
        assert(0 <= l && l <= r && r <= length);
        if (l == r)
            return;
        l += lastnode;
        r += lastnode;
        pushBoundaries(l, r);
        Int a = l, b = r;
        while (a < b) {
            if (a & 1)
                all_apply(a++, f);
            if (b & 1)
                all_apply(--b, f);
            a >>= 1;
            b >>= 1;
        }
        pullBoundaries(l, r);
    }

    // 両端を含む区間にfを作用させる。O(log(n+1))回の演算・作用・合成。
    template <class L, class R> void apply(ClosedSlice<L, R> s, const F &f) {
        apply(resolve_index(length, s.a), resolve_index(length, s.b) + 1, f);
    }

    // f(get(l, r))を満たす最大のrをO(log N)で返します。
    // fは区間の拡大に対して単調で、単位元に対してtrueを返す必要があります。
    template <class Pred> Int max_right(Int l, Pred f) {
        assert(0 <= l && l <= length && f(identity));
        if (l == length)
            return length;
        l += lastnode;
        all_push(l);
        S sm = identity;
        do {
            while (l % 2 == 0)
                l >>= 1;
            if (!f(merge(sm, arr[l]))) {
                while (l < lastnode) {
                    push(l);
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

    // f(get(l, r))を満たす最小のlをO(log N)で返します。
    // fは区間の拡大に対して単調で、単位元に対してtrueを返す必要があります。
    template <class Pred> Int min_left(Int r, Pred f) {
        assert(0 <= r && r <= length && f(identity));
        if (r == 0)
            return 0;
        r += lastnode;
        all_push(r - 1);
        S sm = identity;
        do {
            --r;
            while (r > 1 && r % 2 != 0)
                r >>= 1;
            if (!f(merge(arr[r], sm))) {
                while (r < lastnode) {
                    push(r);
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

// 初期列と演算・作用を指定して構築する。要素数をnとしてO(n+1)回の演算。
template <class S, class F, class V, class Op, class Map, class Comp>
auto initLazySegmentTree(const V &v, Op op, S e, Map map, Comp comp, F id) {
    return LazySegmentTree<S, F>(v, op, e, map, comp, id);
}

// 初期列と演算・作用を指定して構築する。要素数をnとしてO(n+1)回の演算。
template <class V, class Op, class S, class Map, class Comp, class F>
auto newLazySegWith(const V &v, Op op, S e, Map map, Comp comp, F id) {
    return initLazySegmentTree(v, op, e, map, comp, id);
}

#define CPLIB_LAZY_TEMPLATE template <class S, class F, class M, class A, class C, bool B>
#define CPLIB_LAZY_TYPE LazySegmentTree<S, F, M, A, C, B>

// p番目の値を更新する。O(log(n+1))回の演算・作用・合成。
CPLIB_LAZY_TEMPLATE void update(CPLIB_LAZY_TYPE &s, Int p, const S &v) {
    s.update(p, v);
}

// 半開区間[l,r)の積を返す。O(log(n+1))回の演算・作用・合成。
CPLIB_LAZY_TEMPLATE auto get(CPLIB_LAZY_TYPE &s, Int l, Int r) {
    return s.get(l, r);
}

// 全要素の積を返す。O(1)。
CPLIB_LAZY_TEMPLATE auto get_all(const CPLIB_LAZY_TYPE &s) {
    return s.get_all();
}

// 半開区間[l,r)にfを作用させる。O(log(n+1))回の演算・作用・合成。
CPLIB_LAZY_TEMPLATE void apply(CPLIB_LAZY_TYPE &s, Int l, Int r, const F &f) {
    s.apply(l, r, f);
}

// 全要素を空白区切りで文字列化する。O(n log(n+1))回の作用・合成と出力処理。
CPLIB_LAZY_TEMPLATE auto to_string(CPLIB_LAZY_TYPE &s) {
    return s.str();
}

// 両端を含む区間の積を返す。O(log(n+1))回の演算・作用・合成。
template <class S, class F, class M, class A, class C, bool B, class L, class R>
auto get(CPLIB_LAZY_TYPE &s, ClosedSlice<L, R> r) {
    return s.get(r);
}

// 両端を含む区間にfを作用させる。O(log(n+1))回の演算・作用・合成。
template <class S, class F, class M, class A, class C, bool B, class L, class R>
void apply(CPLIB_LAZY_TYPE &s, ClosedSlice<L, R> r, const F &f) {
    s.apply(r, f);
}

// f(積[l,r))が真となる最大のrを返す。単調な述語に対しO(log(n+1))回の演算・作用・合成と述語呼出し。
template <class S, class F, class M, class A, class C, bool B, class P>
auto max_right(CPLIB_LAZY_TYPE &s, Int l, P p) {
    return s.max_right(l, p);
}

// f(積[l,r))が真となる最小のlを返す。単調な述語に対しO(log(n+1))回の演算・作用・合成と述語呼出し。
template <class S, class F, class M, class A, class C, bool B, class P>
auto min_left(CPLIB_LAZY_TYPE &s, Int r, P p) {
    return s.min_left(r, p);
}

#undef CPLIB_LAZY_TEMPLATE
#undef CPLIB_LAZY_TYPE
}
