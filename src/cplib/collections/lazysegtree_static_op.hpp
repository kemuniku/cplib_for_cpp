#pragma once
#include <cplib/collections/private/lazysegtree_base.hpp>

namespace cplib::static_op_lazysegtree {
using detail::lazysegtree_base::apply;
using detail::lazysegtree_base::get;
using detail::lazysegtree_base::get_all;
using detail::lazysegtree_base::max_right;
using detail::lazysegtree_base::min_left;
using detail::lazysegtree_base::to_string;
using detail::lazysegtree_base::update;

// 単位作用を必要としないhasLazy版。演算を静的に特殊化する。
template <class S, class F, class Op, class Map, class Comp>
using LazySegmentTree = detail::lazysegtree_base::LazySegmentTree<S, F, Op, Map, Comp, true>;

// 初期列と静的な演算・作用を指定して構築する。要素数をnとしてO(n+1)回の演算。
template <class S, class F, class V, class Op, class Map, class Comp>
auto initLazySegmentTree(const V &v, Op op, S e, Map map, Comp comp, F id) {
    static_assert(std::is_empty_v<Op> && std::is_empty_v<Map> && std::is_empty_v<Comp>);
    return LazySegmentTree<S, F, Op, Map, Comp>(v, op, e, map, comp, id);
}

// 初期列と静的な演算・作用を指定して構築する。要素数をnとしてO(n+1)回の演算。
template <class V, class Op, class S, class Map, class Comp, class F>
auto newLazySegWith(const V &v, Op op, S e, Map map, Comp comp, F id) {
    return ::cplib::static_op_lazysegtree::initLazySegmentTree(v, op, e, map, comp, id);
}

// 関数として指定した演算・作用と初期列から構築する。要素数をnとしてO(n+1)回の演算。
template <auto Op, auto Map, auto Comp, class V, class S, class F>
auto initLazySegmentTree(const V &v, S e, F id) {
    return ::cplib::static_op_lazysegtree::initLazySegmentTree(
        v, [](S a, S b) { return Op(a, b); }, e, [](F f, S x) { return Map(f, x); },
        [](F f, F g) { return Comp(f, g); }, id);
}
}
