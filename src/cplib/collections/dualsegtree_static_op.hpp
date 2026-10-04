#pragma once
#include <cplib/collections/private/dualsegtree_base.hpp>

namespace cplib {
using detail::dualsegtree_base::apply;
using detail::dualsegtree_base::get;
using detail::dualsegtree_base::len;
using detail::dualsegtree_base::toSeq;
using detail::dualsegtree_base::to_string;
using detail::dualsegtree_base::update;

// 未適用フラグを持ち、作用の単位元を必要としない静的演算版。
template <class S, class F, class Map, class Comp>
using DualSegmentTree = detail::dualsegtree_base::DualSegmentTree<S, F, Map, Comp, true>;

template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(std::span<const S> v, Map map, Comp comp, F id) {
    static_assert(std::is_empty_v<Map> && std::is_empty_v<Comp>);
    return DualSegmentTree<S, F, Map, Comp>(v, map, comp, id);
}

template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(const std::vector<S> &v, Map map, Comp comp, F id) {
    return cplib::initDualSegmentTree(std::span<const S>(v), map, comp, id);
}

template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(Int n, S initial, Map map, Comp comp, F id) {
    static_assert(std::is_empty_v<Map> && std::is_empty_v<Comp>);
    return DualSegmentTree<S, F, Map, Comp>(n, initial, map, comp, id);
}

template <class... Args> auto newDualSegWith(Args &&...args) {
    return cplib::initDualSegmentTree(std::forward<Args>(args)...);
}

template <auto Map, auto Comp, class S, class F>
auto initDualSegmentTree(std::span<const S> v, F id) {
    return cplib::initDualSegmentTree(
        v, [](F f, S x) { return Map(f, x); }, [](F f, F g) { return Comp(f, g); }, id);
}
}
