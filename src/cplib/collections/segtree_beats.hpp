#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <bit>
#include <sstream>

namespace cplib {
// mappingがS.failで失敗を示す場合に子へ降りる汎用Segment Tree Beats。
// 点操作・取得はO(log N)。区間作用の償却計算量はmappingの失敗条件に依存する。
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
    SegmentTreeBeats() = default;

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

    template <class Op, class Map, class Comp>
    SegmentTreeBeats(std::span<const S> values, Op merge, S identity, Map mapping, Comp composition,
                     F id)
        : SegmentTreeBeats(values.size(), merge, identity, mapping, composition, id) {
        std::copy(values.begin(), values.end(), arr.begin() + lastnode_);
        for (Int p = lastnode_ - 1; p; --p)
            pull(p);
    }

    Int len() const {
        return length_;
    }

    Int size() const {
        return len();
    }

    void update(Int p, const S &value) {
        assert(0 <= p && p < length_);
        p += lastnode_;
        all_push(p);
        arr[p] = value;
        for (Int i = 1; i <= log_; ++i)
            pull(p >> i);
    }

    S get(Int p) {
        assert(0 <= p && p < length_);
        all_push(p + lastnode_);
        return arr[p + lastnode_];
    }

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

    template <class L, class R> S get(ClosedSlice<L, R> range) {
        return get(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1);
    }

    template <class L, class R> S operator[](ClosedSlice<L, R> range) {
        return get(range);
    }

    struct Reference {
        SegmentTreeBeats *owner;
        Int index;

        operator S() const {
            return owner->get(index);
        }

        Reference &operator=(const S &x) {
            owner->update(index, x);
            return *this;
        }

        Reference &operator=(const Reference &x) {
            return *this = S(x);
        }
    };

    Reference operator[](Int p) {
        return {this, p};
    }

    Reference operator[](BackwardsIndex p) {
        return (*this)[length_ - p.value];
    }

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

    template <class L, class R> void apply(ClosedSlice<L, R> range, F f) {
        apply(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1, f);
    }

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

template <class S, class F, class V, class Op, class Map, class Comp>
auto initSegmentTreeBeats(const V &v, Op op, S e, Map mapping, Comp composition, F id) {
    return SegmentTreeBeats<S, F>(v, op, e, mapping, composition, id);
}

template <class S, class F, class V, class Op, class Map, class Comp>
void initSegmentTreeBeatsInPlace(SegmentTreeBeats<S, F> &tree, const V &v, Op op, S e, Map mapping,
                                 Comp composition, F id) {
    tree = SegmentTreeBeats<S, F>(v, op, e, mapping, composition, id);
}

template <class S, class F> Int len(const SegmentTreeBeats<S, F> &tree) {
    return tree.len();
}

template <class S, class F> void update(SegmentTreeBeats<S, F> &tree, Int p, const S &value) {
    tree.update(p, value);
}

template <class S, class F> S get(SegmentTreeBeats<S, F> &tree, Int l, Int r) {
    return tree.get(l, r);
}

template <class S, class F, class L, class R>
S get(SegmentTreeBeats<S, F> &tree, ClosedSlice<L, R> range) {
    return tree.get(range);
}

template <class S, class F> void apply(SegmentTreeBeats<S, F> &tree, Int l, Int r, F action) {
    tree.apply(l, r, action);
}

template <class S, class F, class L, class R>
void apply(SegmentTreeBeats<S, F> &tree, ClosedSlice<L, R> range, F action) {
    tree.apply(range, action);
}

template <class S, class F> std::string to_string(SegmentTreeBeats<S, F> &tree) {
    return tree.str();
}

template <class V, class Op, class S, class Map, class Comp, class F>
auto newLazySegWith(const V &v, Op op, S e, Map mapping, Comp composition, F id) {
    return initSegmentTreeBeats(v, op, e, mapping, composition, id);
}
}
