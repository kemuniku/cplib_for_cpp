#pragma once
#include <cplib/collections/private/lazysegtree_base.hpp>

namespace cplib {
template <class T> struct RangeExtremum {
    T value{};
    Int index = -1, left = -1;
    bool operator==(const RangeExtremum &) const = default;
};

template <class T> struct RangeSum {
    T sum{};
    Int len = 0;
    bool operator==(const RangeSum &) const = default;
};

template <class T> struct RangeAffine {
    T a{}, b{};
    bool operator==(const RangeAffine &) const = default;
};

namespace detail {
template <bool Min, bool Assign, class T> auto initRangeExtremumTreeIndex(std::span<const T> v) {
    using S = RangeExtremum<T>;
    auto merge = [](S l, S r) {
        if (l.index < 0)
            return r;
        if (r.index < 0)
            return l;
        S result;
        if constexpr (Min)
            result = r.value < l.value ? r : l;
        else
            result = l.value < r.value ? r : l;
        result.left = l.left;
        return result;
    };
    auto mapping = [](T f, S x) {
        if (x.index < 0)
            return x;
        if constexpr (Assign) {
            x.value = f;
            x.index = x.left;
        } else
            x.value += f;
        return x;
    };
    auto comp = [](T f, T g) {
        if constexpr (Assign)
            return f;
        else
            return f + g;
    };
    std::vector<S> nodes;
    nodes.reserve(v.size());
    for (Int i = 0; i < Int(v.size()); ++i)
        nodes.push_back({v[i], i, i});
    return detail::lazysegtree_base::LazySegmentTree<S, T, decltype(merge), decltype(mapping),
                                                     decltype(comp), true>(
        nodes, merge, S{T{}, -1, -1}, mapping, comp, T{});
}

template <bool Min, bool Assign, class T> auto initRangeExtremumTree(std::span<const T> v) {
    auto merge = [](T l, T r) {
        if constexpr (Min)
            return std::min(l, r);
        else
            return std::max(l, r);
    };
    auto mapping = [](T f, T x) {
        if constexpr (Assign)
            return f;
        else
            return x + f;
    };
    auto comp = [](T f, T g) {
        if constexpr (Assign)
            return f;
        else
            return f + g;
    };
    T identity = Min ? std::numeric_limits<T>::max() : std::numeric_limits<T>::lowest();
    return detail::lazysegtree_base::LazySegmentTree<T, T, decltype(merge), decltype(mapping),
                                                     decltype(comp), true>(v, merge, identity,
                                                                           mapping, comp, T{});
}

template <int Mode, class T> auto initRangeSumTree(std::span<const T> v) {
    using S = RangeSum<T>;
    using F = std::conditional_t<Mode == 2, RangeAffine<T>, T>;
    auto merge = [](S l, S r) { return S{l.sum + r.sum, l.len + r.len}; };
    auto mapping = [](F f, S x) {
        if (x.len == 0)
            return x;
        if constexpr (Mode == 0)
            x.sum += f * T(x.len);
        else if constexpr (Mode == 1)
            x.sum = f * T(x.len);
        else
            x.sum = f.a * x.sum + f.b * T(x.len);
        return x;
    };
    auto comp = [](F f, F g) -> F {
        if constexpr (Mode == 0)
            return f + g;
        else if constexpr (Mode == 1)
            return f;
        else
            return {f.a * g.a, f.a * g.b + f.b};
    };
    std::vector<S> nodes;
    nodes.reserve(v.size());
    for (T x : v)
        nodes.push_back({x, 1});
    return detail::lazysegtree_base::LazySegmentTree<S, F, decltype(merge), decltype(mapping),
                                                     decltype(comp), true>(nodes, merge, S{T{}, 0},
                                                                           mapping, comp, F{});
}
}

// よく使う11種類。構築O(N)、区間更新・取得O(log N)。極値が同じなら最左。
#define CPLIB_RANGE_EXT_FACTORY(NAME, MIN, ASSIGN, INDEX)                                          \
    template <class T> auto NAME(std::span<const T> v) {                                           \
        return detail::initRangeExtremumTree##INDEX<MIN, ASSIGN>(v);                               \
    }                                                                                              \
    template <class T> auto NAME(const std::vector<T> &v) {                                        \
        return NAME(std::span<const T>(v));                                                        \
    }
CPLIB_RANGE_EXT_FACTORY(initRangeAddRangeMinIndex, true, false, Index)
CPLIB_RANGE_EXT_FACTORY(initRangeAddRangeMaxIndex, false, false, Index)
CPLIB_RANGE_EXT_FACTORY(initRangeAssignRangeMinIndex, true, true, Index)
CPLIB_RANGE_EXT_FACTORY(initRangeAssignRangeMaxIndex, false, true, Index)
CPLIB_RANGE_EXT_FACTORY(initRangeAddRangeMin, true, false, )
CPLIB_RANGE_EXT_FACTORY(initRangeAddRangeMax, false, false, )
CPLIB_RANGE_EXT_FACTORY(initRangeAssignRangeMin, true, true, )
CPLIB_RANGE_EXT_FACTORY(initRangeAssignRangeMax, false, true, )
#undef CPLIB_RANGE_EXT_FACTORY
#define CPLIB_RANGE_SUM_FACTORY(NAME, MODE)                                                        \
    template <class T> auto NAME(std::span<const T> v) {                                           \
        return detail::initRangeSumTree<MODE>(v);                                                  \
    }                                                                                              \
    template <class T> auto NAME(const std::vector<T> &v) {                                        \
        return NAME(std::span<const T>(v));                                                        \
    }
CPLIB_RANGE_SUM_FACTORY(initRangeAddRangeSum, 0)
CPLIB_RANGE_SUM_FACTORY(initRangeAssignRangeSum, 1)
CPLIB_RANGE_SUM_FACTORY(initRangeAffineRangeSum, 2)
#undef CPLIB_RANGE_SUM_FACTORY
}
