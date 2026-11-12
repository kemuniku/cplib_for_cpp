#pragma once
#include <cplib/collections/private/segtree_base.hpp>

namespace cplib {
using detail::segtree_base::get;
using detail::segtree_base::get_all;
using detail::segtree_base::len;
using detail::segtree_base::max_right;
using detail::segtree_base::min_left;
using detail::segtree_base::to_string;
using detail::segtree_base::update;

// 演算子の型を静的に保持し、std::functionの間接呼出しを除く独立した特殊化。
// 捕捉なしラムダまたは空の関数オブジェクトを渡す。
template <class T, class Op> using SegmentTree = detail::segtree_base::SegmentTree<T, Op>;

// 長さnの単位元列から木を構築する。O(n+1)回の演算。
template <class T, class Op> auto initSegmentTree(Int n, Op op, T e) {
    static_assert(std::is_empty_v<Op>, "static op requires a stateless functor");
    return SegmentTree<T, Op>(n, op, e);
}

// 配列から木を構築する。要素数をnとしてO(n+1)回の演算。
template <class T, class Op> auto initSegmentTree(std::span<const T> v, Op op, T e) {
    static_assert(std::is_empty_v<Op>, "static op requires a stateless functor");
    return SegmentTree<T, Op>(v, op, e);
}

// 配列から木を構築する。要素数をnとしてO(n+1)回の演算。
template <class T, class Op> auto initSegmentTree(const std::vector<T> &v, Op op, T e) {
    return cplib::initSegmentTree(std::span<const T>(v), op, e);
}

// 配列から木を構築する。要素数をnとしてO(n+1)回の演算。
template <class V, class Op, class T> auto newSegWith(const V &v, Op op, T e) {
    return cplib::initSegmentTree(v, op, e);
}

// 名前付き関数は非型テンプレート引数で指定できる。
// 静的に指定した演算と配列から木を構築する。要素数をnとしてO(n+1)回の演算。
template <auto Op, class V, class T> auto initSegmentTree(const V &v, T e) {
    return cplib::initSegmentTree(v, [](T a, T b) { return Op(a, b); }, e);
}
}
