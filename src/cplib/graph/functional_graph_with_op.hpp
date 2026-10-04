#pragma once
#include <cplib/graph/functional_graph.hpp>
#include <cplib/collections/segtree.hpp>

namespace cplib {
// HLDの逆順とサイクル順を別々のセグメント木に保持し、非可換な積の順序を保つ。
template <class T> class FunctionalGraph_with_op : public FunctionalGraph {
    std::function<T(T, T)> op_;
    T e_;
    SegmentTree<T> st_hld_, st_cycle_;
    std::vector<Int> offsets_;

    Int index(Int v) const {
        return tree.N - 1 - tree.toSeq(v);
    }

    T value(Int v) const {
        return st_hld_[index(v)];
    }

public:
    template <class Op>
    FunctionalGraph_with_op(FunctionalGraph graph, const std::vector<T> &values, Op op, T e)
        : FunctionalGraph(std::move(graph)), op_(op), e_(e), st_hld_(0, op, e),
          st_cycle_(0, op, e) {
        assert(values.size() == cycle_number.size());
        std::vector<T> hld(tree.N, e), cyc;
        for (Int v = 0; v < Int(values.size()); ++v)
            hld[index(v)] = values[v];
        st_hld_ = SegmentTree<T>(hld, op, e);
        for (const auto &vertices : cycle) {
            offsets_.push_back(cyc.size());
            for (Int v : vertices)
                cyc.push_back(values[v]);
        }
        st_cycle_ = SegmentTree<T>(cyc, op, e);
    }

    void set(Int x, const T &value) {
        st_hld_.update(index(x), value);
        if (incycle(x))
            st_cycle_.update(offsets_[cycle_number[x]] + cycle_idx[x], value);
    }

    T get(Int x) const {
        return value(x);
    }

    struct Reference {
        FunctionalGraph_with_op *owner;
        Int index;

        operator T() const {
            return owner->get(index);
        }

        Reference &operator=(const T &v) {
            owner->set(index, v);
            return *this;
        }

        Reference &operator=(const Reference &v) {
            return *this = T(v);
        }
    };

    Reference operator[](Int x) {
        return {this, x};
    }

    T operator[](Int x) const {
        return value(x);
    }

    // 可換・冪等な積について全始点の到達可能頂点の積をO(N)で求める。
    std::vector<T> prod_reachable_idempotent_all() const {
        Int n = cycle_number.size();
        std::vector<T> result(n, e_);
        for (const auto &vertices : cycle) {
            T product = e_;
            for (Int v : vertices)
                product = op_(product, value(v));
            for (Int v : vertices)
                result[v] = product;
        }
        for (Int v : tree.I)
            if (v < n && !incycle(v))
                result[v] = op_(value(v), result[tree.P[v]]);
        return result;
    }

    // k回の移動で訪れる値の積。木のHLD分解と周回積の二乗法でO(log² N+log k)。
    // startからk回移動するまでの積。include_start=falseなら始点を積に含めない。
    T prod(Int start, Int k, bool include_start = true) const {
        assert(k >= 0);
        if (!include_start) {
            if (!k)
                return e_;
            return prod(movekth(start, 1), k - 1);
        }
        T result = e_;
        Int r = roots[start], last = r;
        bool ends_in_tree = depth(start) > k;
        if (ends_in_tree)
            last = movekth(start, k);
        for (auto [l, u] : tree.path(last, start, ends_in_tree, true))
            result = op_(result, st_hld_.get(l, u));
        if (ends_in_tree)
            return result;
        Int cid = cycle_number[r], size = cyclesize(start), offset = offsets_[cid],
            begin = cycle_idx[r];
        UInt count = UInt(k - depth(start)) + 1, full = count / size, remainder = count % size;
        T product = op_(st_cycle_.get(offset + begin, offset + size),
                        st_cycle_.get(offset, offset + begin));
        for (; full; full >>= 1) {
            if (full & 1)
                result = op_(result, product);
            if (full > 1)
                product = op_(product, product);
        }
        Int first = std::min<Int>(remainder, size - begin);
        result = op_(result, st_cycle_.get(offset + begin, offset + begin + first));
        if (UInt(first) < remainder)
            result = op_(result, st_cycle_.get(offset, offset + remainder - first));
        return result;
    }

    // 累積積の列は最初の積を計算した後に1頂点ずつ延長する。O(log² N+log l+r-l)。
    // prod(start,l), ..., prod(start,r-1)を配列で返す。O(log^2 N + log l + (r-l))。
    std::vector<T> prod_range(Int start, Int l, Int r, bool include_start = true) const {
        assert(0 <= l && l <= r);
        if (l == r)
            return {};
        std::vector<T> result;
        result.reserve(r - l);
        result.push_back(prod(start, l, include_start));
        if (r - l == 1)
            return result;
        Int now = movekth(start, l + 1);
        for (Int i = l + 1; i < r; ++i) {
            result.push_back(op_(result.back(), value(now)));
            now = next(now);
        }
        return result;
    }

    // prod(start,l), ..., prod(start,r-1)を順にfoldで畳み込む。O(log^2 N + log l + (r-l))。
    template <class Fold, class U>
    U prod_range_fold(Int start, Int l, Int r, Fold fold, U e, bool include_start = true) const {
        assert(0 <= l && l <= r);
        if (l == r)
            return e;
        T prefix = prod(start, l, include_start);
        U result = fold(e, prefix);
        if (r - l == 1)
            return result;
        Int now = movekth(start, l + 1);
        for (Int i = l + 1; i < r; ++i) {
            prefix = op_(prefix, value(now));
            result = fold(result, prefix);
            now = next(now);
        }
        return result;
    }

    // 単調な判定の最初の失敗位置（移動回数）を求める。O(log² N+log L)。
    // xから始め、累積積に対するpredicateが初めてfalseになるまでの移動距離を返す。
    // 移動距離の上限はlimitとし、limit回移動してもtrueならlimitを返す。
    template <class Predicate> Int move_while(Predicate predicate, Int x, Int limit) const {
        assert(0 <= limit && limit < std::numeric_limits<Int>::max());
        Int count = limit + 1, used = 1;
        T accumulated = value(x);
        if (!predicate(accumulated))
            return 0;
        if (used == count)
            return limit;
        Int r = roots[x];
        bool first_segment = true;
        for (auto [l, u] : tree.path(r, x, true, true)) {
            Int begin = l + Int(first_segment);
            first_segment = false;
            Int end = begin + std::min(u - begin, count - used);
            if (begin < end) {
                Int bound =
                    st_hld_.max_right(begin, [&](T v) { return predicate(op_(accumulated, v)); });
                if (bound < end)
                    return used + bound - begin;
                accumulated = op_(accumulated, st_hld_.get(begin, end));
                used += end - begin;
            }
            if (used == count)
                return limit;
        }
        Int cid = cycle_number[r], size = cyclesize(x), offset = offsets_[cid],
            begin = (cycle_idx[r] + 1) % size;
        auto range = [&](Int l, Int u) { return st_cycle_.get(offset + l, offset + u); };
        T one_cycle = op_(range(begin, size), range(0, begin));
        Int max_cycles = (count - used) / size;
        if (max_cycles) {
            std::vector<T> powers{one_cycle};
            for (Int block = 1; block <= max_cycles / 2; block *= 2)
                powers.push_back(op_(powers.back(), powers.back()));
            Int accepted = 0;
            for (Int i = Int(powers.size()) - 1; i >= 0; --i) {
                Int cycles = Int(1) << i;
                if (cycles <= max_cycles - accepted) {
                    T candidate = op_(accumulated, powers[i]);
                    if (predicate(candidate)) {
                        accumulated = candidate;
                        accepted += cycles;
                    }
                }
            }
            used += accepted * size;
            if (used == count)
                return limit;
        }
        auto consume = [&](Int l, Int u) {
            if (l == u)
                return true;
            Int bound = st_cycle_.max_right(offset + l,
                                            [&](T v) { return predicate(op_(accumulated, v)); });
            if (bound < offset + u) {
                used += bound - offset - l;
                return false;
            }
            accumulated = op_(accumulated, range(l, u));
            used += u - l;
            return true;
        };
        Int rest = std::min(count - used, size), first = std::min(rest, size - begin);
        if (!consume(begin, begin + first))
            return used;
        rest -= first;
        if (rest && !consume(0, rest))
            return used;
        return used - 1;
    }
};

template <class Source, class T, class Op>
auto initFunctionalGraph_with_op(const Source &source, const std::vector<T> &values, Op op, T e) {
    if constexpr (std::is_same_v<Source, FunctionalGraph>)
        return FunctionalGraph_with_op<T>(source, values, op, e);
    else
        return FunctionalGraph_with_op<T>(initFunctionalGraph(source), values, op, e);
}

template <class T> void set(FunctionalGraph_with_op<T> &f, Int x, const T &value) {
    f.set(x, value);
}

template <class T> auto prod_reachable_idempotent_all(const FunctionalGraph_with_op<T> &f) {
    return f.prod_reachable_idempotent_all();
}

template <class T>
T prod(const FunctionalGraph_with_op<T> &f, Int x, Int k, bool include_start = true) {
    return f.prod(x, k, include_start);
}

template <class T>
auto prod_range(const FunctionalGraph_with_op<T> &f, Int x, Int l, Int r,
                bool include_start = true) {
    return f.prod_range(x, l, r, include_start);
}

template <class T, class F, class U>
U prod_range_fold(const FunctionalGraph_with_op<T> &graph, Int x, Int l, Int r, F f, U e,
                  bool include_start = true) {
    return graph.prod_range_fold(x, l, r, f, e, include_start);
}

template <class T, class F>
Int move_while(const FunctionalGraph_with_op<T> &graph, F f, Int x, Int limit) {
    return graph.move_while(f, x, limit);
}
}
