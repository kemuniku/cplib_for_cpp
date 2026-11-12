#pragma once
#include <cplib/collections/private/range_reverse_treap.hpp>

namespace cplib {
template <class S, class F>
using RangeReverseDualSegmentTree = detail::RangeReverseTreap<S, F, false, true>;

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class S, class F, class Map, class Comp>
auto initRangeReverseDualSegmentTree(std::span<const S> values, Map mapping, Comp composition,
                                     F id) {
    return RangeReverseDualSegmentTree<S, F>(values, {}, S{}, mapping, composition, id);
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class S, class F, class Map, class Comp>
auto initRangeReverseDualSegmentTree(const std::vector<S> &values, Map mapping, Comp composition,
                                     F id) {
    return RangeReverseDualSegmentTree<S, F>(values, {}, S{}, mapping, composition, id);
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class S, class F, class Map, class Comp>
auto initRangeReverseDualSegmentTree(Int n, S value, Map mapping, Comp composition, F id) {
    assert(n >= 0);
    return initRangeReverseDualSegmentTree(std::vector<S>(n, value), mapping, composition, id);
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class V, class Map, class Comp, class F>
auto newRangeReverseDualSegWith(const V &values, Map mapping, Comp composition, F id) {
    return initRangeReverseDualSegmentTree(values, mapping, composition, id);
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class S, class Map, class Comp, class F>
auto newRangeReverseDualSegWith(Int n, S value, Map mapping, Comp composition, F id) {
    return initRangeReverseDualSegmentTree(n, value, mapping, composition, id);
}
}
