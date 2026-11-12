#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <bit>
#include <sstream>

namespace cplib {
// mappingがS.failで失敗を示す場合に子へ降りる汎用Segment Tree Beats。
// 各演算は基本処理O(log(N+1))に、mapping失敗で子へ降りた回数を加えた計算量。
// 償却計算量はmappingの失敗条件に依存する。各コールバックはO(1)とする。
template <class S, class F> class SegmentTreeBeats {
    S identity_{};
    F id_{};
    std::function<S(S, S)> merge_;
    std::function<S(F, S)> mapping_;
    std::function<F(F, F)> composition_;
    Int lastnode_ = 1, log_ = 0, length_ = 0;

    void pull(Int p) {
        arr[p] = merge_(arr[2 * p], arr[2 * p + 1]);
    }

    void all_apply(Int p, F f) {
        arr[p] = mapping_(f, arr[p]);
        if (p < lastnode_) {
            lazy[p] = composition_(f, lazy[p]);
            if (arr[p].fail) {
                push(p);
                pull(p);
            }
        }
    }

    void push(Int p) {
        all_apply(2 * p, lazy[p]);
        all_apply(2 * p + 1, lazy[p]);
        lazy[p] = id_;
    }

    void all_push(Int p) {
        for (Int i = log_; i >= 1; --i)
            push(p >> i);
    }

    void push_boundaries(Int l, Int r) {
        for (Int i = log_; i > std::countr_zero(UInt(l)); --i)
            push(l >> i);
        for (Int i = log_; i > std::countr_zero(UInt(r)); --i)
            push((r - 1) >> i);
    }

public:
    std::vector<S> arr;
    std::vector<F> lazy;
    // 未初期化の木を作る。O(1)。
    SegmentTreeBeats() = default;

    // n個の単位元から木を構築する。O(n+1)。
    template <class Op, class Map, class Comp>
    SegmentTreeBeats(Int n, Op merge, S identity, Map mapping, Comp composition, F id)
        : identity_(identity), id_(id), merge_(merge), mapping_(mapping), composition_(composition),
          length_(n) {
        assert(n >= 0);
        while (lastnode_ < n) {
            lastnode_ *= 2;
            ++log_;
        }
        arr.assign(2 * lastnode_, identity);
        lazy.assign(lastnode_, id);
    }

    // 配列から木を構築する。演算がO(1)ならO(N+1)。
    template <class Op, class Map, class Comp>
    SegmentTreeBeats(std::span<const S> values, Op merge, S identity, Map mapping, Comp composition,
                     F id)
        : SegmentTreeBeats(values.size(), merge, identity, mapping, composition, id) {
        std::copy(values.begin(), values.end(), arr.begin() + lastnode_);
        for (Int p = lastnode_ - 1; p; --p)
            pull(p);
    }

    // 要素数を返す。O(1)。
    Int len() const {
        return length_;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return len();
    }

    // 位置pの値を更新する。O(log(N+1)+D)、Dは遅延伝播時に失敗して子へ降りた回数。
    void update(Int p, const S &value) {
        assert(0 <= p && p < length_);
        p += lastnode_;
        all_push(p);
        arr[p] = value;
        for (Int i = 1; i <= log_; ++i)
            pull(p >> i);
    }

    // 位置pの値を返す。O(log(N+1)+D)、Dは遅延伝播時に失敗して子へ降りた回数。
    S get(Int p) {
        assert(0 <= p && p < length_);
        all_push(p + lastnode_);
        return arr[p + lastnode_];
    }

    // O(log(N+1)+D)、Dは遅延伝播時に失敗して子へ降りた回数。
    // 半開区間[l,r)の演算結果を返します。
    S get(Int l, Int r) {
        assert(0 <= l && l <= r && r <= length_);
        if (l == r)
            return identity_;
        l += lastnode_;
        r += lastnode_;
        push_boundaries(l, r);
        S left = identity_, right = identity_;
        while (l < r) {
            if (l & 1)
                left = merge_(left, arr[l++]);
            if (r & 1)
                right = merge_(arr[--r], right);
            l >>= 1;
            r >>= 1;
        }
        return merge_(left, right);
    }

    // 閉区間を集約する。O(log(N+1)+D)、Dは遅延伝播時に失敗して子へ降りた回数。
    template <class L, class R> S get(ClosedSlice<L, R> range) {
        return get(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1);
    }

    // 閉区間を集約する。O(log(N+1)+D)、Dは遅延伝播時に失敗して子へ降りた回数。
    template <class L, class R> S operator[](ClosedSlice<L, R> range) {
        return get(range);
    }

    struct Reference {
        SegmentTreeBeats *owner;
        Int index;

        // 現在値を読み取る。O(log(N+1)+D)。
        operator S() const {
            return owner->get(index);
        }

        // 値を代入する。O(log(N+1)+D)。
        Reference &operator=(const S &x) {
            owner->update(index, x);
            return *this;
        }

        // 別の参照の値を読み取り、この位置へ代入する。O(log(N+1)+D)。
        Reference &operator=(const Reference &x) {
            return *this = S(x);
        }
    };

    // 位置pへの読み書き用参照を返す。取得O(1)、読み書きの計算量はget/updateと同じ。
    Reference operator[](Int p) {
        return {this, p};
    }

    // 末尾から指定した位置への読み書き用参照を返す。取得O(1)、読み書きの計算量はget/updateと同じ。
    Reference operator[](BackwardsIndex p) {
        return (*this)[length_ - p.value];
    }

    // O(log(N+1)+D)、Dはmappingが失敗して子へ降りた回数。
    // 半開区間[l,r)にfを作用させます。
    void apply(Int l, Int r, F f) {
        assert(0 <= l && l <= r && r <= length_);
        if (l == r)
            return;
        l += lastnode_;
        r += lastnode_;
        push_boundaries(l, r);
        for (Int a = l, b = r; a < b; a >>= 1, b >>= 1) {
            if (a & 1)
                all_apply(a++, f);
            if (b & 1)
                all_apply(--b, f);
        }
        for (Int i = std::countr_zero(UInt(l)) + 1; i <= log_; ++i)
            pull(l >> i);
        for (Int i = std::countr_zero(UInt(r)) + 1; i <= log_; ++i)
            pull((r - 1) >> i);
    }

    // 閉区間にfを作用させる。O(log(N+1)+D)、Dはmappingが失敗して子へ降りた回数。
    template <class L, class R> void apply(ClosedSlice<L, R> range, F f) {
        apply(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1, f);
    }

    // 全要素を文字列で返す。出力長L、各取得の失敗回数の合計DとしてO(N log(N+1)+D+L)。
    std::string str() {
        std::ostringstream out;
        for (Int i = 0; i < length_; ++i) {
            if (i)
                out << ' ';
            out << get(i);
        }
        return out.str();
    }
};

// 初期配列または要素数から木を構築する。演算がO(1)ならO(N+1)。
template <class S, class F, class V, class Op, class Map, class Comp>
auto initSegmentTreeBeats(const V &v, Op op, S e, Map mapping, Comp composition, F id) {
    return SegmentTreeBeats<S, F>(v, op, e, mapping, composition, id);
}

// 既存の木を再構築する。元の要素数をM、新しい要素数をNとし演算がO(1)ならO(M+N+1)。
template <class S, class F, class V, class Op, class Map, class Comp>
void initSegmentTreeBeatsInPlace(SegmentTreeBeats<S, F> &tree, const V &v, Op op, S e, Map mapping,
                                 Comp composition, F id) {
    tree = SegmentTreeBeats<S, F>(v, op, e, mapping, composition, id);
}

// 要素数を返す。O(1)。
template <class S, class F> Int len(const SegmentTreeBeats<S, F> &tree) {
    return tree.len();
}

// 点を更新する。O(log(N+1)+D)、Dは伝播中のmapping失敗回数。
template <class S, class F> void update(SegmentTreeBeats<S, F> &tree, Int p, const S &value) {
    tree.update(p, value);
}

// 半開区間を集約する。O(log(N+1)+D)、Dは伝播中のmapping失敗回数。
template <class S, class F> S get(SegmentTreeBeats<S, F> &tree, Int l, Int r) {
    return tree.get(l, r);
}

// 閉区間を集約する。O(log(N+1)+D)、Dは伝播中のmapping失敗回数。
template <class S, class F, class L, class R>
S get(SegmentTreeBeats<S, F> &tree, ClosedSlice<L, R> range) {
    return tree.get(range);
}

// 半開区間へ作用させる。O(log(N+1)+D)、Dはmapping失敗回数。
template <class S, class F> void apply(SegmentTreeBeats<S, F> &tree, Int l, Int r, F action) {
    tree.apply(l, r, action);
}

// 閉区間へ作用させる。O(log(N+1)+D)、Dはmapping失敗回数。
template <class S, class F, class L, class R>
void apply(SegmentTreeBeats<S, F> &tree, ClosedSlice<L, R> range, F action) {
    tree.apply(range, action);
}

// 全要素を文字列で返す。出力長L、取得時のmapping失敗回数の合計DとしてO(N log(N+1)+D+L)。
template <class S, class F> std::string to_string(SegmentTreeBeats<S, F> &tree) {
    return tree.str();
}

// 初期配列または要素数から木を構築する。演算がO(1)ならO(N+1)。
template <class V, class Op, class S, class Map, class Comp, class F>
auto newLazySegWith(const V &v, Op op, S e, Map mapping, Comp composition, F id) {
    return initSegmentTreeBeats(v, op, e, mapping, composition, id);
}
}
