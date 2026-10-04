#pragma once
#include <cplib/graph/functional_graph.hpp>
#include <cplib/collections/lazysegtree.hpp>

namespace cplib {
// サイクル頂点をcycle側だけに保持する、遅延作用付きの独立した実装。
// STはvalue_type/action_type、calc_e/op/id/compositionと遅延木の操作を提供する。
template <class ST> class FunctionalGraph_with_lazy_op : public FunctionalGraph {
public:
    using S = typename ST::value_type;
    using F = typename ST::action_type;

private:
    ST st_hld_, st_cycle_;
    std::vector<Int> offsets_;

    Int hld_index(Int v) const {
        return tree.N - 1 - tree.toSeq(v);
    }

    Int cycle_index(Int v) const {
        return offsets_[cycle_number[v]] + cycle_idx[v];
    }

    S op(const S &a, const S &b) const {
        return st_hld_.calc_op(a, b);
    }

    S identity() const {
        return st_hld_.calc_e();
    }

    F action_power(F base, Int count) const {
        F result = st_hld_.calc_id();
        for (; count; count >>= 1) {
            if (count & 1)
                result = st_hld_.calc_composition(base, result);
            if (count > 1)
                base = st_hld_.calc_composition(base, base);
        }
        return result;
    }

    void apply_cycle_prefix(Int cid, Int begin, Int count, const F &action) {
        if (!count)
            return;
        Int size = cycle[cid].size(), offset = offsets_[cid], first = std::min(count, size - begin);
        st_cycle_.apply(offset + begin, offset + begin + first, action);
        if (first < count)
            st_cycle_.apply(offset, offset + count - first, action);
    }

    std::vector<S> walk_values(Int start, Int count) {
        std::vector<S> out;
        out.reserve(count);
        if (!count)
            return out;
        Int tree_count = std::min(count, depth(start));
        if (tree_count) {
            Int last = movekth(start, tree_count - 1);
            for (auto [l, r] : tree.path(last, start, true, true))
                st_hld_.append_range(l, r, out);
        }
        Int cycle_count = count - tree_count;
        if (!cycle_count)
            return out;
        Int r = roots[start], cid = cycle_number[r], size = cycle[cid].size(),
            offset = offsets_[cid], begin = cycle_idx[r];
        std::vector<S> one_cycle;
        Int materialize_count = std::min(cycle_count, size),
            first = std::min(materialize_count, size - begin);
        st_cycle_.append_range(offset + begin, offset + begin + first, one_cycle);
        if (first < materialize_count)
            st_cycle_.append_range(offset, offset + materialize_count - first, one_cycle);
        for (Int i = 0; i < cycle_count; ++i)
            out.push_back(one_cycle[i % one_cycle.size()]);
        return out;
    }

public:
    template <class Factory>
    FunctionalGraph_with_lazy_op(FunctionalGraph graph, const std::vector<S> &values, Factory make)
        : FunctionalGraph(std::move(graph)), st_hld_(make(std::vector<S>{})),
          st_cycle_(make(std::vector<S>{})) {
        assert(values.size() == cycle_number.size());
        std::vector<S> hld(tree.N, identity()), cyc;
        for (Int v = 0; v < Int(values.size()); ++v)
            if (!incycle(v))
                hld[hld_index(v)] = values[v];
        for (const auto &vertices : cycle) {
            offsets_.push_back(cyc.size());
            for (Int v : vertices)
                cyc.push_back(values[v]);
        }
        st_hld_ = make(hld);
        st_cycle_ = make(cyc);
    }

    // 頂点xの現在値を返す。O(log N)
    S get(Int x) {
        return incycle(x) ? st_cycle_.get(cycle_index(x)) : st_hld_.get(hld_index(x));
    }

    // 頂点xの値をvalueに変更する。O(log N)
    void set(Int x, const S &value) {
        if (incycle(x))
            st_cycle_.update(cycle_index(x), value);
        else
            st_hld_.update(hld_index(x), value);
    }

    struct Reference {
        FunctionalGraph_with_lazy_op *owner;
        Int index;

        operator S() const {
            return owner->get(index);
        }

        Reference &operator=(const S &v) {
            owner->set(index, v);
            return *this;
        }

        Reference &operator=(const Reference &v) {
            return *this = S(v);
        }
    };

    Reference operator[](Int x) {
        return {this, x};
    }

    // 同じ頂点への複数訪問は作用の合成の二乗法で処理する。O(log² N+log k)。
    // startからk回移動するまでに訪れる各頂点へactionを作用させる。
    // 同じ頂点を複数回訪れた場合は、その回数だけactionを作用させる。
    // include_start=falseなら始点には作用させない。O(log^2 N + log k)。
    void apply(Int start, Int k, const F &action, bool include_start = true) {
        assert(0 <= k && k < std::numeric_limits<Int>::max());
        if (!include_start) {
            if (k)
                apply(movekth(start, 1), k - 1, action);
            return;
        }
        Int d = depth(start);
        if (k < d) {
            for (auto [l, r] : tree.path(movekth(start, k), start, true, true))
                st_hld_.apply(l, r, action);
            return;
        }
        Int r = roots[start];
        for (auto [l, u] : tree.path(r, start, false, true))
            st_hld_.apply(l, u, action);
        Int cid = cycle_number[r], size = cycle[cid].size(), count = k - d + 1, full = count / size,
            remainder = count % size, offset = offsets_[cid];
        if (full)
            st_cycle_.apply(offset, offset + size, action_power(action, full));
        apply_cycle_prefix(cid, cycle_idx[r], remainder, action);
    }

    // 先に全遅延値を伝播し、可換・冪等な積を全頂点についてO(N)で求める。
    std::vector<S> prod_reachable_idempotent_all() {
        Int n = cycle_number.size();
        std::vector<S> result(n, identity());
        auto hld = st_hld_.materialize(), cyc = st_cycle_.materialize();
        auto value = [&](Int v) { return incycle(v) ? cyc[cycle_index(v)] : hld[hld_index(v)]; };
        for (const auto &vertices : cycle) {
            S product = identity();
            for (Int v : vertices)
                product = op(product, value(v));
            for (Int v : vertices)
                result[v] = product;
        }
        for (Int v : tree.I)
            if (v < n && !incycle(v))
                result[v] = op(value(v), result[tree.P[v]]);
        return result;
    }

    // startからk回移動するまでの積。include_start=falseなら始点を積に含めない。
    // O(log^2 N + log k)
    S prod(Int start, Int k, bool include_start = true) {
        assert(0 <= k && k < std::numeric_limits<Int>::max());
        if (!include_start) {
            if (!k)
                return identity();
            return prod(movekth(start, 1), k - 1);
        }
        S result = identity();
        Int r = roots[start], last = r;
        bool ends_in_tree = depth(start) > k;
        if (ends_in_tree)
            last = movekth(start, k);
        for (auto [l, u] : tree.path(last, start, ends_in_tree, true))
            result = op(result, st_hld_.get(l, u));
        if (ends_in_tree)
            return result;
        Int cid = cycle_number[r], size = cycle[cid].size(), offset = offsets_[cid],
            begin = cycle_idx[r], count = k - depth(start) + 1, full = count / size,
            remainder = count % size;
        if (full) {
            S product = op(st_cycle_.get(offset + begin, offset + size),
                           st_cycle_.get(offset, offset + begin));
            for (; full; full >>= 1) {
                if (full & 1)
                    result = op(result, product);
                if (full > 1)
                    product = op(product, product);
            }
        }
        Int first = std::min(remainder, size - begin);
        if (first)
            result = op(result, st_cycle_.get(offset + begin, offset + begin + first));
        if (first < remainder)
            result = op(result, st_cycle_.get(offset, offset + remainder - first));
        return result;
    }

    // 遅延木の区間を一括展開するため、O(log² N+log l+r-l)を保つ。
    // prod(start,l), ..., prod(start,r-1)を配列で返す。
    // O(log^2 N + log l + (r-l))。
    std::vector<S> prod_range(Int start, Int l, Int r, bool include_start = true) {
        assert(0 <= l && l <= r);
        if (l == r)
            return {};
        std::vector<S> result;
        result.reserve(r - l);
        result.push_back(prod(start, l, include_start));
        if (r - l == 1)
            return result;
        auto values = walk_values(movekth(start, l + 1), r - l - 1);
        for (const S &value : values)
            result.push_back(op(result.back(), value));
        return result;
    }

    // prod(start,l), ..., prod(start,r-1)を順にfoldで畳み込む。
    // O(log^2 N + log l + (r-l))。
    template <class Fold, class U>
    U prod_range_fold(Int start, Int l, Int r, Fold fold, U e, bool include_start = true) {
        assert(0 <= l && l <= r);
        if (l == r)
            return e;
        S prefix = prod(start, l, include_start);
        U result = fold(e, prefix);
        if (r - l == 1)
            return result;
        auto values = walk_values(movekth(start, l + 1), r - l - 1);
        for (const S &value : values) {
            prefix = op(prefix, value);
            result = fold(result, prefix);
        }
        return result;
    }

    // xから始め、prefixのprodに対するpredicateが初めてfalseになる移動距離を返す。
    // limit回移動してもtrueならlimitを返す。O(log^2 N + log limit)。
    // predicate(e)=trueであり、一度falseになった後は頂点を追加してもfalseのまま、
    // というACL max_rightと同じ単調性を仮定する。
    template <class Predicate> Int move_while(Predicate predicate, Int x, Int limit) {
        assert(0 <= limit && limit < std::numeric_limits<Int>::max());
        Int count = limit + 1, used = 1;
        S accumulated = get(x);
        if (!predicate(accumulated))
            return 0;
        if (used == count)
            return limit;
        Int r = roots[x];
        bool first_segment = true;
        for (auto [l, u] : tree.path(r, x, false, true)) {
            Int begin = l + Int(first_segment);
            first_segment = false;
            Int end = begin + std::min(u - begin, count - used);
            if (begin < end) {
                Int bound =
                    st_hld_.max_right(begin, [&](S v) { return predicate(op(accumulated, v)); });
                if (bound < end)
                    return used + bound - begin;
                accumulated = op(accumulated, st_hld_.get(begin, end));
                used += end - begin;
            }
            if (used == count)
                return limit;
        }
        Int cid = cycle_number[r], size = cycle[cid].size(), offset = offsets_[cid],
            begin = incycle(x) ? (cycle_idx[x] + 1) % size : cycle_idx[r];
        auto range = [&](Int l, Int u) { return st_cycle_.get(offset + l, offset + u); };
        S one_cycle = op(range(begin, size), range(0, begin));
        Int max_cycles = (count - used) / size;
        if (max_cycles) {
            std::vector<S> powers{one_cycle};
            for (Int block = 1; block <= max_cycles / 2; block *= 2)
                powers.push_back(op(powers.back(), powers.back()));
            Int accepted = 0;
            for (Int i = Int(powers.size()) - 1; i >= 0; --i) {
                Int cycles = Int(1) << i;
                if (cycles <= max_cycles - accepted) {
                    S candidate = op(accumulated, powers[i]);
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
            Int bound =
                st_cycle_.max_right(offset + l, [&](S v) { return predicate(op(accumulated, v)); });
            if (bound < offset + u) {
                used += bound - offset - l;
                return false;
            }
            accumulated = op(accumulated, range(l, u));
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

// 関数オブジェクト版。e/idはACLと同じゼロ引数関数。
template <class Source, class S, class Op, class E, class Map, class Comp, class Id>
auto initFunctionalGraph_with_lazy_op(const Source &source, const std::vector<S> &values, Op op,
                                      E e, Map mapping, Comp composition, Id id) {
    using F = std::invoke_result_t<Id>;
    using ST = LazySegmentTree<S, F>;
    auto make = [&](const std::vector<S> &v) { return ST(v, op, e(), mapping, composition, id()); };
    if constexpr (std::is_same_v<Source, FunctionalGraph>)
        return FunctionalGraph_with_lazy_op<ST>(source, values, make);
    else
        return FunctionalGraph_with_lazy_op<ST>(initFunctionalGraph(source), values, make);
}

// 明示した型はvectorから構築でき、上記のSTインターフェイスを提供する。
template <class ST, class Source>
auto initFunctionalGraph_with_lazy_op(const Source &source,
                                      const std::vector<typename ST::value_type> &values) {
    auto make = [](const auto &v) { return ST(v); };
    if constexpr (std::is_same_v<Source, FunctionalGraph>)
        return FunctionalGraph_with_lazy_op<ST>(source, values, make);
    else
        return FunctionalGraph_with_lazy_op<ST>(initFunctionalGraph(source), values, make);
}

template <class Source, class S, class Op, class E, class Map, class Comp, class Id>
auto initFunctionalGraph_with_op(const Source &source, const std::vector<S> &values, Op op, E e,
                                 Map mapping, Comp composition, Id id) {
    return initFunctionalGraph_with_lazy_op(source, values, op, e, mapping, composition, id);
}

template <class ST> auto get(FunctionalGraph_with_lazy_op<ST> &graph, Int x) {
    return graph.get(x);
}

template <class ST>
void set(FunctionalGraph_with_lazy_op<ST> &graph, Int x, const typename ST::value_type &v) {
    graph.set(x, v);
}

template <class ST> auto prod_reachable_idempotent_all(FunctionalGraph_with_lazy_op<ST> &graph) {
    return graph.prod_reachable_idempotent_all();
}

template <class ST>
auto prod(FunctionalGraph_with_lazy_op<ST> &graph, Int x, Int k, bool include_start = true) {
    return graph.prod(x, k, include_start);
}

template <class ST>
void apply(FunctionalGraph_with_lazy_op<ST> &graph, Int x, Int k,
           const typename ST::action_type &action, bool include_start = true) {
    graph.apply(x, k, action, include_start);
}

template <class ST>
auto prod_range(FunctionalGraph_with_lazy_op<ST> &graph, Int x, Int l, Int r,
                bool include_start = true) {
    return graph.prod_range(x, l, r, include_start);
}

template <class ST, class Fold, class U>
auto prod_range_fold(FunctionalGraph_with_lazy_op<ST> &graph, Int x, Int l, Int r, Fold f, U e,
                     bool include_start = true) {
    return graph.prod_range_fold(x, l, r, f, e, include_start);
}

template <class ST, class Predicate>
Int move_while(FunctionalGraph_with_lazy_op<ST> &graph, Predicate f, Int x, Int limit) {
    return graph.move_while(f, x, limit);
}
}
