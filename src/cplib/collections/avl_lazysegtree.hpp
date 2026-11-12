#pragma once
#include <cplib/collections/private/avl_sequence.hpp>

namespace cplib {
template <class S, class F> using AvlLazySegmentTree = detail::AvlSequence<S, F, true, true>;
// 集約・区間作用付きAVL列を作る。mappingは集約値にも作用すること。時間・領域 O(N)。
template <class S, class F, class Op, class Map, class Comp>
auto initAvlLazySegmentTree(const std::vector<S> &values, Op op, S e, Map mapping, Comp composition, F id) {
    return AvlLazySegmentTree<S, F>(values, op, e, mapping, composition, id);
}
// 単位元をn個持つ集約・区間作用付きAVL列を作る。時間・領域 O(n)。
template <class S, class F, class Op, class Map, class Comp>
auto initAvlLazySegmentTree(Int n, Op op, S e, Map mapping, Comp composition, F id) {
    assert(n >= 0);
    return initAvlLazySegmentTree(std::vector<S>(n, e), op, e, mapping, composition, id);
}
// 配列または長さから集約・区間作用付きAVL列を作る。時間・領域 O(N)。
template <class V, class Op, class S, class Map, class Comp, class F>
auto newAvlLazySegWith(const V &values, Op op, S e, Map mapping, Comp composition, F id) {
    return initAvlLazySegmentTree(values, op, e, mapping, composition, id);
}
}
