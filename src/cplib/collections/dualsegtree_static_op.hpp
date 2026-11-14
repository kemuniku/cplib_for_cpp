#pragma once
#include <cplib/collections/private/dualsegtree_base.hpp>

namespace cplib::static_op_dualsegtree {
using detail::dualsegtree_base::apply;
using detail::dualsegtree_base::get;
using detail::dualsegtree_base::toSeq;
using detail::dualsegtree_base::to_string;
using detail::dualsegtree_base::update;

// 未適用フラグを持ち、作用の単位元を必要としない静的演算版。
template <class S, class F, class Map, class Comp>
using DualSegmentTree = detail::dualsegtree_base::DualSegmentTree<S, F, Map, Comp, true>;

// 配列と静的な作用・合成を指定して構築する。要素数をnとしてO(n+1)。
template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(std::span<const S> v, Map map, Comp comp, F id) {
    static_assert(std::is_empty_v<Map> && std::is_empty_v<Comp>);
    return DualSegmentTree<S, F, Map, Comp>(v, map, comp, id);
}

// 配列と静的な作用・合成を指定して構築する。要素数をnとしてO(n+1)。
template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(const std::vector<S> &v, Map map, Comp comp, F id) {
    return ::cplib::static_op_dualsegtree::initDualSegmentTree(std::span<const S>(v), map, comp, id);
}

// 長さnの同じ初期値を持つ列から構築する。O(n+1)。
template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(Int n, S initial, Map map, Comp comp, F id) {
    static_assert(std::is_empty_v<Map> && std::is_empty_v<Comp>);
    return DualSegmentTree<S, F, Map, Comp>(n, initial, map, comp, id);
}

// 初期列と作用・合成を指定して構築する。要素数をnとしてO(n+1)。
template <class... Args> auto newDualSegWith(Args &&...args) {
    return ::cplib::static_op_dualsegtree::initDualSegmentTree(std::forward<Args>(args)...);
}

// 配列と関数として指定した作用・合成から構築する。要素数をnとしてO(n+1)。
template <auto Map, auto Comp, class S, class F>
auto initDualSegmentTree(std::span<const S> v, F id) {
    return ::cplib::static_op_dualsegtree::initDualSegmentTree(
        v, [](F f, S x) { return Map(f, x); }, [](F f, F g) { return Comp(f, g); }, id);
}
}
